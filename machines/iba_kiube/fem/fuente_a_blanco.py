"""Simulación de punta a punta (v0.35): fuente PIG → región central → aceleración → lámina extractora → blanco.

Todas las etapas ven el MISMO campo magnético: el mapa FE oficial (out/fem_mag/ngsolve) × coil.escala_iso, con
f_RF = 40.640 MHz (el nivel isócrono de v0.34).

  1. rc     Región central con el menisco autoconsistente del plasma (fem/central_region_menisco.py, puller con bisel,
            chimenea Ø6) y el campo E de la malla FE; las partículas salen a r = 126 mm (≈1.3 MeV) o tras `--vueltas`.
            El estado final (posición, momento por ajuste cuadrático de los 3 últimos puntos, energía, fase RF) se
            entrega a la etapa 2.
  2. acel   fem/aceleracion.acelerar desde ese estado hasta la lámina (azimut del stripper del puerto, r ≥ r de la
            órbita cerrada de 18 MeV en ese azimut, 10 mm de lámina).
  3. ext    Cambio de carga H⁻ → p (misma velocidad), dispersión múltiple en 5 µm de carbono (Highland), trayectoria
            de los protones en el mapa hasta r = 0.60 m y recta hasta la ventana del blanco (cara del yugo + 10 mm).
            El azimut de la lámina se reajusta (propagando el H⁻ hacia atrás o adelante, sin RF en ese tramo) para
            centrar el haz en el blanco del puerto, como hace el operador con el carrusel.
Uso: python fem/fuente_a_blanco.py rc|acel|ext|todo [--n 800] [--puerto 0]  → out/fuente_blanco/
"""
from __future__ import annotations
import argparse, json, math, pathlib, sys, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
OUT = ROOT / "out/fuente_blanco"
NPZ = ROOT / "out/fem_mag/ngsolve/bz_midplane_NI115k.npz"
F_RF = 40.640e6
C = 299792458.0; QE = 1.602176634e-19
M_H_MEV = 939.294; M_P_MEV = 938.272
MEN = ROOT / "out/central_region/menisco_ch5/menisco_jp1440.json"


SUF = ""; RC = "rc_estados.npz"
ESCALA = None                                                        # --escala: nivel de campo distinto de coil.escala_iso


def params():
    from kiube import load_params
    P, _ = load_params()
    return P, float(ESCALA or P["coil"].get("escala_iso", 1.0))


# ----------------------------------------------------------------------------------------------- 1. región central
def etapa_rc(n=800, vueltas=16.0, semilla=7):
    import central_region as cr
    P, k = params()
    M = json.load(open(MEN))
    niv = sorted(M["niveles"], key=lambda r: -r["cos"])
    cs = np.array([r["cos"] for r in niv]); bordes = [(cs[i] + cs[i + 1]) / 2 for i in range(len(cs) - 1)]
    pr0 = dict(cr.PRESETS["kiube"]); pr0["_nombre"] = "kiube"; pr0["puller_slit_h"] = 10.0
    pr0["menisco_s_cav"] = pr0["chimney_r"] - 1.2
    pr0.update(B_source="fem_map", B_npz=str(NPZ), B_escala=k, f_rf=F_RF)
    rng = np.random.default_rng(semilla)
    PH = rng.uniform(-80.0, 40.0, n)                      # ventana de emisión del caso base (fuera de ella no hay haz útil)
    k_of = np.array([int(np.sum(c < np.array(bordes))) for c in np.cos(np.radians(PH))])
    fil = dict(PH=[], motivo=[], X=[], V=[], T=[], t=[], nivel=[])
    for kk, r in enumerate(niv):
        sel = np.flatnonzero(k_of == kk)
        if not len(sel):
            continue
        pr = dict(pr0, _menisco=r["coef"], _PH=PH[sel].tolist())
        mesh, gfu, info, st = cr.solve(pr, None)
        Ef = cr.FEField(mesh, gfu, False)
        res, trs, B, Bf = cr.run_species(Ef, pr, info["fuentes"][0], P, (-70.0, 30.0, 10.0), vueltas, len(sel),
                                         0.9 * pr["R_dom"])
        print(f"[rc] cos {r['cos']}: {len(sel)} → {res['perdidas_haz']}  ({Bf.src})", flush=True)
        for b in B:
            tr = b["tr"]
            fil["PH"].append(b["fase"]); fil["motivo"].append(b["motivo"]); fil["nivel"].append(r["cos"])
            if len(tr) >= 3:
                t3 = tr[-3:, 0]; x3 = tr[-3:, 1:4]
                # derivada en el último punto del polinomio de Lagrange por los 3 últimos puntos (paso no uniforme)
                t0_, t1_, t2_ = t3
                d0 = (t2_ - t1_) / ((t0_ - t1_) * (t0_ - t2_)); d1 = (t2_ - t0_) / ((t1_ - t0_) * (t1_ - t2_))
                d2 = (2 * t2_ - t0_ - t1_) / ((t2_ - t0_) * (t2_ - t1_))
                v = (d0 * x3[0] + d1 * x3[1] + d2 * x3[2]) * 1e-3 / 1e-9          # m/s
                fil["X"].append((tr[-1, 1:4] * 1e-3).tolist()); fil["V"].append(v.tolist())
                fil["T"].append(float(tr[-1, 4]) / 1e3); fil["t"].append(float(tr[-1, 0]))
            else:
                fil["X"].append([0, 0, 0]); fil["V"].append([0, 0, 0]); fil["T"].append(0.0); fil["t"].append(0.0)
        del mesh, gfu, Ef
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / RC, **{k_: np.array(v) for k_, v in fil.items()}, escala=k, n=n)
    m = np.array(fil["motivo"])
    print("[rc]", {x: int((m == x).sum()) for x in set(m)}, flush=True)


def _estado_acel(d):
    """Estados de la etapa 1 que siguen vivos o salieron por r → entrada de acelerar (momento con |p| de T)."""
    m = d["motivo"]; ok = np.isin(m, ["vive", "r"]) & (d["T"] > 0.2)
    X = d["X"][ok]; V = d["V"][ok]; T = d["T"][ok]
    g = 1 + T / M_H_MEV; p = g * M_H_MEV * 1e6 * QE / C ** 2 * C * np.sqrt(1 - 1 / g ** 2)
    Pv = V / np.linalg.norm(V, axis=1)[:, None] * p[:, None]
    w = 2 * math.pi * F_RF
    phi0 = np.radians(d["PH"][ok]) + w * d["t"][ok] * 1e-9            # cos(w t' + φ0) con t' = 0 al empezar la etapa 2
    return ok, dict(X=X, P=Pv, T=T, phi0=phi0)


def lamina(M, P, puerto):
    """Azimut del stripper del puerto (óptica de v0.2, centrado en el puerto medido) y radio interior = órbita cerrada
    de 18 MeV en ese azimut en el mapa FE."""
    import orbitas as O
    from kiube.beam import extraction
    E = extraction(P)
    s = E["strippers"][puerto]
    phi = s["phi"] % 360
    thm = (phi - 45.0) % 90.0 + 45.0
    o = O.orbita_cerrada(M, 18.0, 0.45, th0=math.radians(thm))
    return dict(phi_deg=phi, r_in=o["r0_m"], r_out=o["r0_m"] + P["extraction"]["foil_size"] / 1000), E


# ----------------------------------------------------------------------------------------------- 2. aceleración
def etapa_acel(puerto=0, max_vueltas=400):
    import orbitas as O, aceleracion as A
    P, k = params()
    d = dict(np.load(OUT / RC, allow_pickle=True))
    ok, est = _estado_acel(d)
    M = O.Mapa(str(NPZ), escala=k)
    lam, _ = lamina(M, P, puerto)
    print(f"[acel] {ok.sum()} partículas desde la región central; lámina {lam}", flush=True)
    t0 = time.time()
    R = A.acelerar(M, None, d["PH"][ok], 0, 0, f_rf=F_RF, estado=est, lamina=lam, max_vueltas=max_vueltas,
                   log=lambda s: print(s, flush=True))
    R.update(lamina=lam, escala=k, puerto=puerto, t_s=time.time() - t0, idx_rc=np.flatnonzero(ok).tolist())
    json.dump(R, open(OUT / f"acel_p{puerto}{SUF}.json", "w"), default=float)
    print("[acel]", R["motivos"], f"{time.time() - t0:.0f} s", flush=True)


# ----------------------------------------------------------------------------------------------- 3. extracción
class MapaExt:
    """Campo para los protones extraídos: mapa principal (r < 0.575 m, armónicos) + mapa extendido de la FE
    (bz_midplane_ext: 0.40–0.95 m, con el material de cada punto) por interpolación bilineal periódica de 90°.
    Puntos en hierro (yugo) → B = 0 y la partícula se da por perdida en el yugo (en_hierro)."""
    R_UNION = 0.575

    def __init__(self, M, escala):
        from scipy.interpolate import RegularGridInterpolator as RGI
        d = np.load(NPZ.parent / NPZ.name.replace("bz_midplane_", "bz_midplane_ext_"))
        self.M = M; self.k = escala
        r = d["r_m"]; ph = d["phi_deg"]; B = d["Bz_T"]
        self.fe = np.isnan(B) | (d["material"] == "iron")
        Bz = np.where(self.fe, 0.0, np.nan_to_num(B))
        self.ib = RGI((r, ph), Bz, bounds_error=False, fill_value=0.0)
        self.ih = RGI((r, ph), self.fe.astype(float), method="nearest", bounds_error=False, fill_value=0.0)
        self.th0 = M.th0

    def _ph(self, th):
        return (np.degrees(th) - 45.0) % 90.0 + 45.0

    def B(self, r, th):
        r = np.asarray(r, float); th = np.asarray(th, float)
        Bz, dr, dt = self.M.B(np.minimum(r, self.R_UNION), th)
        ext = r >= self.R_UNION
        if ext.any():
            re, pe = r[ext], self._ph(th[ext])
            f = lambda a, b: self.ib(np.c_[a, b]) * self.k
            Bz[ext] = f(re, pe)
            dr[ext] = (f(re + 0.002, pe) - f(re - 0.002, pe)) / 0.004
            dt[ext] = (f(re, pe + 0.5) - f(re, pe - 0.5)) / np.radians(1.0)
        return Bz, dr, dt

    def en_hierro(self, X):
        r = np.hypot(X[:, 0], X[:, 1]); th = np.arctan2(X[:, 1], X[:, 0])
        return (r > 0.6) & (self.ih(np.c_[r, self._ph(th)]) > 0.5)


def _boris(M, X, Pm, q, m_kg, dt, pasos=None, parar=None):
    import aceleracion as A
    X = X.copy(); Pm = Pm.copy()
    vivo = np.ones(len(X), bool); k = 0
    while vivo.any() and (pasos is None or k < pasos):
        i = np.flatnonzero(vivo)
        B = A.B3(M, X[i])
        gm = np.sqrt(1 + np.einsum("ij,ij->i", Pm[i], Pm[i]) / (m_kg * C) ** 2)
        tv = (q * dt / (2 * m_kg * gm))[:, None] * B
        pp = Pm[i] + np.cross(Pm[i], tv)
        sv = 2 * tv / (1 + np.einsum("ij,ij->i", tv, tv))[:, None]
        Pm[i] = Pm[i] + np.cross(pp, sv)
        gm = np.sqrt(1 + np.einsum("ij,ij->i", Pm[i], Pm[i]) / (m_kg * C) ** 2)
        X[i] += Pm[i] / (m_kg * gm)[:, None] * dt
        k += 1
        if parar is not None:
            vivo[i[parar(X[i])]] = False
    return X, Pm


def _a_lamina(M, X, Pm, dphi_deg):
    """Propaga el H⁻ (sin RF) un ángulo azimutal dphi (negativo = hacia atrás en el tiempo)."""
    if abs(dphi_deg) < 1e-9:
        return X, Pm
    m = M_H_MEV * 1e6 * QE / C ** 2
    th0 = np.degrees(np.arctan2(X[:, 1], X[:, 0]))
    dt = 0.02e-9 * np.sign(dphi_deg)
    obj = th0 + dphi_deg
    Xo = X.copy(); Po = Pm.copy(); hecho = np.zeros(len(X), bool)
    Xc, Pc = X.copy(), Pm.copy()
    for _ in range(200000):
        if hecho.all():
            break
        i = np.flatnonzero(~hecho)
        Xn, Pn = _boris(M, Xc[i], Pc[i], -QE, m, dt, pasos=1)
        th = np.degrees(np.arctan2(Xn[:, 1], Xn[:, 0]))
        d = (th - obj[i] + 180) % 360 - 180
        cruzo = (d >= 0) if dphi_deg > 0 else (d <= 0)
        Xc[i], Pc[i] = Xn, Pn
        j = i[cruzo]; Xo[j] = Xn[cruzo]; Po[j] = Pn[cruzo]; hecho[j] = True
    return Xo, Po


def protones_a_blanco(M, Xs, Ps, P, E, puerto, semilla=11):
    """H⁻ en la lámina → protones; trayectoria hasta r = 0.60 m y recta hasta el plano de la ventana del blanco."""
    rng = np.random.default_rng(semilla)
    mp = M_P_MEV * 1e6 * QE / C ** 2
    Pp = Ps * (M_P_MEV / M_H_MEV)                                     # misma velocidad
    pmag = np.linalg.norm(Pp, axis=1)
    g = np.sqrt(1 + (pmag / (mp * C)) ** 2); beta = np.sqrt(1 - 1 / g ** 2)
    pc = pmag * C / QE / 1e6
    xX0 = P["extraction"]["stripper_foil_um"] * 1e-4 * 2.0 / 42.7    # carbono 2.0 g/cm³, X0 = 42.7 g/cm²
    th0 = 13.6 / (beta * pc) * math.sqrt(xX0) * (1 + 0.038 * math.log(xX0))
    u = Pp / pmag[:, None]
    e1 = np.cross(u, [0, 0, 1.0]); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(e1, u)
    a1 = rng.normal(0, th0, len(u)); a2 = rng.normal(0, th0, len(u))
    u = u + a1[:, None] * e1 + a2[:, None] * e2; u /= np.linalg.norm(u, axis=1)[:, None]
    Pp = u * pmag[:, None]
    ME = MapaExt(M, M.escala)
    H = P["yoke"]["side"] / 2000
    hierro = np.zeros(len(Xs), bool)

    def parar(X_):
        fuera = (np.abs(X_[:, 0]) > H) | (np.abs(X_[:, 1]) > H)
        return fuera | ME.en_hierro(X_)
    X, Pf = _boris(ME, Xs, Pp, QE, mp, 0.05e-9, parar=parar)
    hierro = ME.en_hierro(X) & ~((np.abs(X[:, 0]) > H) | (np.abs(X[:, 1]) > H))
    # recta hasta el plano del blanco del puerto: centro y dirección de diseño (cara del yugo + 10 mm)
    t = E["targets"][puerto]
    c = np.array([t["x"], t["y"], 0.0]) * 1e-3; n_ = np.array([math.cos(math.radians(t["dir"])), math.sin(math.radians(t["dir"])), 0])
    d = Pf / np.linalg.norm(Pf, axis=1)[:, None]
    s = ((c - X) @ n_) / (d @ n_)
    Y = X + s[:, None] * d
    eh = np.cross([0, 0, 1.0], n_)                                    # horizontal en el plano del blanco
    xh = (Y - c) @ eh; zv = Y[:, 2]
    ang_h = np.arcsin(np.clip(d @ eh, -1, 1)) * 1e3; ang_v = np.arcsin(np.clip(d[:, 2], -1, 1)) * 1e3
    Tp = (g - 1) * M_P_MEV
    return dict(x_mm=xh * 1e3, z_mm=zv * 1e3, xp_mrad=ang_h, zp_mrad=ang_v, T_MeV=Tp, theta0_mrad=th0 * 1e3, yugo=hierro)


def etapa_ext(puerto=0, dphis=np.arange(-3.0, 3.01, 0.5)):
    import orbitas as O
    P, k = params()
    R = json.load(open(OUT / f"acel_p{puerto}{SUF}.json"))
    M = O.Mapa(str(NPZ), escala=k)
    lam, E = lamina(M, P, puerto)
    est = R["estados_lamina"]
    ids = sorted(est, key=int)
    Xs = np.array([est[i]["X"] for i in ids]); Ps = np.array([est[i]["P"] for i in ids])
    Ts = np.array([est[i]["T"] for i in ids]); fs = np.array([est[i]["fase"] for i in ids])
    print(f"[ext] {len(ids)} H⁻ en la lámina: T {Ts.mean():.3f} ± {Ts.std():.3f} MeV", flush=True)
    filas = []
    for dp in dphis:
        X1, P1 = _a_lamina(M, Xs, Ps, dp)
        r = protones_a_blanco(M, X1, P1, P, E, puerto)
        filas.append(dict(dphi=float(dp), x_medio=float(r["x_mm"].mean()), xp_medio=float(r["xp_mrad"].mean())))
        print(f"[ext] Δφ lámina {dp:+.1f}°: centro del haz en el blanco {r['x_mm'].mean():+.1f} mm, "
              f"ángulo {r['xp_mrad'].mean():+.1f} mrad", flush=True)
    xs = np.array([f["x_medio"] for f in filas]); dps = np.array([f["dphi"] for f in filas])
    o = np.argsort(xs); dbest = float(np.interp(0.0, xs[o], dps[o])) if xs.min() < 0 < xs.max() else float(dps[np.argmin(np.abs(xs))])
    X1, P1 = _a_lamina(M, Xs, Ps, dbest)
    r = protones_a_blanco(M, X1, P1, P, E, puerto)
    col = P["extraction"]["collimator_d"] / 2; win = P["target"]["window_d"] / 2
    rr = np.where(r["yugo"], np.inf, np.hypot(r["x_mm"], r["z_mm"]))
    res = dict(puerto=puerto, lamina=lam, dphi_lamina_ajuste_deg=dbest, barrido=filas, n=len(ids),
               T_MeV=[float(r["T_MeV"].mean()), float(r["T_MeV"].std())],
               sigma_x_mm=float(r["x_mm"].std()), sigma_z_mm=float(r["z_mm"].std()),
               sigma_xp_mrad=float(r["xp_mrad"].std()), sigma_zp_mrad=float(r["zp_mrad"].std()),
               centro_mm=[float(r["x_mm"].mean()), float(r["z_mm"].mean())],
               frac_colimador=float(np.mean(rr <= col)), frac_ventana=float(np.mean(rr <= win)),
               frac_yugo=float(np.mean(r["yugo"])),
               theta0_lamina_mrad=float(np.mean(r["theta0_mrad"])), fase_lamina=[float(fs.mean()), float(fs.std())],
               puntos=dict(x=r["x_mm"].tolist(), z=r["z_mm"].tolist(), T=r["T_MeV"].tolist(),
                           xp=r["xp_mrad"].tolist(), zp=r["zp_mrad"].tolist()))
    json.dump(res, open(OUT / f"ext_p{puerto}{SUF}.json", "w"), default=float, indent=1)
    print("[ext]", json.dumps({k_: v for k_, v in res.items() if k_ not in ("puntos", "barrido")}, default=float), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("etapa", choices=["rc", "acel", "ext", "todo"])
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--vueltas", type=float, default=16.0)
    ap.add_argument("--puerto", type=int, default=0)
    ap.add_argument("--escala", type=float, default=None)
    ap.add_argument("--suf", default="")
    ap.add_argument("--rc", default="rc_estados.npz")
    a = ap.parse_args()
    ESCALA = a.escala; SUF = a.suf; RC = a.rc
    if a.etapa in ("rc", "todo"):
        etapa_rc(a.n, a.vueltas)
    if a.etapa in ("acel", "todo"):
        etapa_acel(a.puerto)
    if a.etapa in ("ext", "todo"):
        etapa_ext(a.puerto)
