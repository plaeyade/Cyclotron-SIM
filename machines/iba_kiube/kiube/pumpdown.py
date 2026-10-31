"""Pump-down y presión de operación del KIUBE con TPMC (física de Molflow+) sobre el STL del vacío.

Fuentes de gas (todas citadas en params.yaml):
  - H₂O de metales sin hornear: q(t) = q1h · (1 h / t)^α · f_material      [CERN-CH Ec.13, Tabla 3]
  - H₂O de juntas de Viton: desgasificación del elastómero expuesto q10h·(10 h/t)^0.5  [CERN-CH Fig.30]
                            + permeación del agua del aire, que crece hasta el régimen estacionario [CERN-CH Sec.6.2]
  - H₂ de la fuente PIG activa (operación)
Respuesta lineal: una corrida TPMC por grupo de material da R_i = P / q_i (mbar por mbar·l/s/cm²);
P(t) = Σ_i q_i(t) R_i + Σ_juntas Q_j(t) R_j + P0  (cuasi-estático: V/S ≈ 0.1 s ≪ t).
"""
from __future__ import annotations

import json, math, pathlib
import numpy as np


def _load(outdir):
    from .tpmc import load_gas_mesh
    d = pathlib.Path(outdir)
    G = json.load(open(d / "kiube_molflow_grupos.json"))["grupos"]
    m = load_gas_mesh(d / "kiube_vacio_cm.stl")
    return m, G


def _mask(G, n, names):
    msk = np.zeros(n, bool)
    for g in names:
        a, b = G[g]["facetas"]
        msk[a - 1:b] = True
    return msk


def seal_permeation_fraction(t_h, t_ss_days=60.0):
    """Fracción del flujo estacionario de permeación a tiempo t (lámina, Crank): 1 + 2Σ(−1)^n e^(−n²t/τ)."""
    tau = t_ss_days * 24.0 / 5.3
    t = np.atleast_1d(np.asarray(t_h, float))
    s = np.ones_like(t)
    for n in range(1, 60):
        s += 2 * (-1) ** n * np.exp(-n * n * t / tau)
    return np.clip(s, 0, 1)


def run(P, outdir, n_particles=30000, seed=11):
    from .tpmc import TPMC, sticking_for_speed
    V = P["vacuum"]
    m, G = _load(outdir)
    nF = len(m.faces)
    c = m.triangles_center * 10.0                                  # cm → mm
    A = m.area_faces
    r = np.hypot(c[:, 0], c[:, 1])
    beam = (r < 500) & (np.abs(c[:, 2]) < 70) & _mask(G, nF, [g for g in G if g.startswith("pared_")])
    pumps = [g for g in G if g.startswith("bomba_") and g in V["active_pumps"]]
    out = {"grupos_fuente": {}, "zona_haz": "facetas de pared con r < 500 mm y |z| < 70 mm (polos, valles, dees)"}

    def sticking(mass, S):
        s = np.zeros(nF)
        for g in pumps:
            a, b = G[g]["facetas"]
            s[a - 1:b] = sticking_for_speed(S, G[g]["area_cm2"], mass)
        return s

    def mean_beam(Pf):
        return float((Pf * A)[beam].sum() / A[beam].sum())

    def mean_pumps(Pf):
        mk = _mask(G, nF, pumps)
        return float((Pf * A)[mk].sum() / A[mk].sum())

    # --- H2O: una corrida por material (q = 1 mbar·l/s/cm² → respuesta R)
    sH2O = sticking(18.015, V["pump_speed_H2O_l_s"])
    t = TPMC(m, sH2O, seed=seed)
    resp = {}
    for g in sorted(G):
        if not g.startswith(("pared_", "lamina_")) or G[g]["area_cm2"] < 1:
            continue
        w = np.zeros(nF); a, b = G[g]["facetas"]; w[a - 1:b] = A[a - 1:b]
        Q = float(w.sum())                                          # caudal con q = 1
        res = t.run(w, n_particles=n_particles)
        Pf = t.pressure(res, Q, 18.015)
        resp[g] = dict(area_cm2=Q, R_haz=mean_beam(Pf), R_bombas=mean_pumps(Pf),
                       S_eff_l_s=Q / mean_beam(Pf), rebotes=res["mean_bounces"], perdidas=res["lost"])
        out["grupos_fuente"][g] = resp[g]
    # juntas: se reparten sobre la pared de inox de la cámara (anillo exterior) → respuesta de pared_inox
    R_seal = resp["pared_inox"]["R_haz"] / resp["pared_inox"]["area_cm2"]     # mbar por mbar·l/s

    # --- H2 de la fuente activa (operación)
    sH2 = sticking(2.016, V["pump_speed_l_s"])
    t2 = TPMC(m, sH2, seed=seed + 1)
    gsrc = f"fuente_h2_{V['active_source']}"
    w = np.zeros(nF); a, b = G[gsrc]["facetas"]; w[a - 1:b] = A[a - 1:b]
    Qh2 = V["h2_sccm_per_active_source"] * 0.01689                 # sccm → mbar·l/s
    res = t2.run(w, n_particles=n_particles)
    Pf2 = t2.pressure(res, Qh2, 2.016)
    out["H2_fuente"] = dict(Q_mbar_l_s=Qh2, P_haz_mbar=mean_beam(Pf2), P_bombas_mbar=mean_pumps(Pf2),
                            S_eff_l_s=Qh2 / mean_beam(Pf2), rebotes=res["mean_bounces"])
    # perfil radial de H2 en el plano medio (para pérdidas de H⁻)
    edges = np.arange(0, 560, 20.0)
    prof = []
    for r0, r1 in zip(edges[:-1], edges[1:]):
        mk = beam & (r >= r0) & (r < r1)
        prof.append([0.5 * (r0 + r1), float((Pf2 * A)[mk].sum() / A[mk].sum()) if mk.any() else float("nan")])
    out["H2_perfil_radial"] = prof

    # --- pump-down
    th = np.logspace(math.log10(0.5), math.log10(2000), 120)
    q1, al = V["outgassing_q1h_H2O"], V["outgassing_alpha"]
    f = V["outgassing_factor"]
    P_met = sum(f.get(g, 1.0) * q1 * th ** (-al) * resp[g]["R_haz"] for g in resp)
    L_cm = V["seal_length_m"] * 100.0
    Q_perm = V["seal_perm_mbar_l_s_cm"] * L_cm * seal_permeation_fraction(th)
    Q_elast = V["polymer_q10h_H2O"] * (10.0 / th) ** 0.5 * V["seal_exposed_cm2"]
    P_seal = (Q_perm + Q_elast) * R_seal
    P_tot = P_met + P_seal + V["P0_ultimate_mbar"]
    out["pumpdown"] = dict(t_h=th.tolist(), P_metales=P_met.tolist(), P_juntas=P_seal.tolist(), P_total=P_tot.tolist())
    Pm = V["P_base_measured_mbar"]
    reach = th[np.argmax(P_tot <= Pm)] if (P_tot <= Pm).any() else None
    out["comparacion"] = dict(P_base_medida_mbar=Pm, t_para_alcanzarla_h=None if reach is None else float(reach),
                              P_10h=float(np.interp(10, th, P_tot)), P_24h=float(np.interp(24, th, P_tot)),
                              P_100h=float(np.interp(100, th, P_tot)), P_1000h=float(np.interp(1000, th, P_tot)),
                              S_eff_H2O_l_s=float(sum(resp[g]["area_cm2"] for g in resp) / sum(resp[g]["R_haz"] for g in resp)))
    return out


def curve(P, T, f_acero=1.0, t_h=None):
    """P(t) [mbar] en la zona del haz a partir de las respuestas TPMC guardadas."""
    V = P["vacuum"]
    th = np.asarray(T["pumpdown"]["t_h"] if t_h is None else t_h, float)
    q1, al, f = V["outgassing_q1h_H2O"], V["outgassing_alpha"], dict(V["outgassing_factor"])
    f["pared_acero"] = f_acero * f.get("pared_acero", 1.0)
    R = T["grupos_fuente"]
    P_met = sum(f.get(g, 1.0) * q1 * th ** (-al) * R[g]["R_haz"] for g in R)
    R_seal = R["pared_inox"]["R_haz"] / R["pared_inox"]["area_cm2"]
    Q_perm = V["seal_perm_mbar_l_s_cm"] * V["seal_length_m"] * 100.0 * seal_permeation_fraction(th)
    Q_el = V["polymer_q10h_H2O"] * (10.0 / th) ** 0.5 * V["seal_exposed_cm2"]
    return th, P_met + (Q_perm + Q_el) * R_seal + V["P0_ultimate_mbar"]


def plot(P, T, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    V = P["vacuum"]
    th = np.asarray(T["pumpdown"]["t_h"])
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    fac = V["outgassing_factor"].get("pared_acero", 1.0)
    for fa, ls, lab in ((1, "-", f"calibrado: acero del yugo ×{fac:g} (bombeo TUD03T)"),
                        (1 / fac, ":", "acero del yugo = inox (CERN-CH, metal limpio)")):
        _, Pt = curve(P, T, fa)
        ax.loglog(th, Pt, "k" + ls if ls != ":" else "k:", lw=1.8 if fa == 1 else 1.2, label=lab)
    pm = V.get("pumpdown_medido")
    if pm:
        ax.loglog(np.array(pm["t_min"]) / 60, pm["P_mbar"], "rv", ms=7, label="medido: 4 × VHS-6 tras 2 h abierto (TUD03T)")
    ax.loglog(th, T["pumpdown"]["P_metales"], color="#3b6fd8", lw=1, label="solo metales (calibrado)")
    ax.loglog(th, T["pumpdown"]["P_juntas"], color="#e05a2b", lw=1, label="juntas de Viton (desgasificación + permeación)")
    ax.axhline(V["P_base_measured_mbar"], color="#c0392b", ls="-.", lw=1.2, label="KIUBE base medida, con RF (TUD03T: 2.7e-7)")
    ax.axhline(5.3e-7, color="#8e44ad", ls="-.", lw=1.0, label="Cyclone KEY base medida (KEY22)")
    ax.axhline(T["H2_fuente"]["P_haz_mbar"], color="#16a085", ls="--", lw=1.0,
               label=f"H₂ con fuente encendida ({V['h2_sccm_per_active_source']:g} sccm), TPMC")
    ax.set_xlabel("tiempo de bombeo [h]"); ax.set_ylabel("presión en la zona del haz [mbar]")
    ax.set_title("KIUBE: pump-down (TPMC sobre el STL de la cámara; acero calibrado con TUD03T)")
    ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=7.5, loc="lower left")
    ax.set_ylim(1e-10, 1e-4)
    fig.tight_layout(); fig.savefig(path, dpi=140)
