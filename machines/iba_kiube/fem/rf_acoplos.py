"""Acoplador de potencia (lazo inductivo) del resonador del KIUBE, calculado con el modo propio de fem/rf_kiube.py.

El lazo (plano que contiene el eje del stem, junto a la tapa = máximo de corriente del resonador λ/4) recoge el flujo
magnético del modo: fem = ω·|Φ|. Con el modo normalizado a V_dee = 40 kV y el lazo cargado por la línea de Z0 = 50 Ω,
la potencia que el resonador entregaría a la línea es P_ext = fem² / (2·Z0) y Q_ext = ω·W / P_ext.
Acoplamiento crítico (β = 1, sin potencia reflejada con el haz a plena corriente):
    P_ext(40 kV) = P_paredes + P_haz,   P_paredes = P_cobre / factor_Q_real,   P_haz = T·I_haz
La reactancia propia del lazo (ωL ≈ 10–30 Ω) se compensa con la longitud de la línea o un condensador serie: no cambia
el acoplamiento a primer orden.
"""
from __future__ import annotations
import math
import numpy as np
import ngsolve as ng

MU0 = 4e-7 * math.pi


def flujo_lazo(g, s, f_MHz, P, h_max_mm=200.0, n=48):
    """Densidad de flujo normal al lazo en una rejilla (ρ, z) bajo la tapa de la carcasa del dee 0.
    Devuelve (rho_m, z_m desde la tapa hacia abajo, Bn [T] a 40 kV) con el modo g escalado por s."""
    mesh = g.space.mesh
    om = 2 * math.pi * f_MHz * 1e6
    B = ng.curl(g) * (1j / om)                                           # B = rot E / (−iω)
    DE, Y, X, RF = P["dee"], P["yoke"], P["exterior"], P["rf"]
    a = math.radians(DE["valleys"][0])
    hx, hy = Y["pump_hole_r_pos"] * math.cos(a) / 1e3, Y["pump_hole_r_pos"] * math.sin(a) / 1e3
    al = a + math.radians(RF["acoplador_az"]); ux, uy = math.cos(al), math.sin(al)
    nx, ny = -uy, ux                                                     # normal al plano del lazo (azimutal)
    zcap = (Y["height"] / 2 + X["rf_housing_h"] - 20) / 1e3
    r_in = DE["stem_radius"] / 1e3; r_out = (X["rf_housing_d"] / 2 - 20) / 1e3
    rho = np.linspace(r_in + 2e-3, r_out - 2e-3, n)
    dz = np.linspace(1e-3, h_max_mm / 1e3, n)
    Bn = np.zeros((n, n))
    for i, r in enumerate(rho):
        for j, d in enumerate(dz):
            try:
                v = B(mesh(hx + r * ux, hy + r * uy, zcap - d))
                Bn[i, j] = abs(v[0] * nx + v[1] * ny) * s
            except Exception:                                            # noqa: BLE001
                Bn[i, j] = np.nan
    return rho, dz, Bn


def acoplador(res, g, s, P, log=print):
    """Altura del lazo [mm] para acoplamiento crítico con la extensión radial rf.acoplador_rho."""
    m = next(mm for mm in res["modos"] if mm["fase"] == "en fase")
    f, W, Pcu = m["f_MHz"], m["W_J_40kV"], m["P_kW_40kV"] * 1e3
    RF = P["rf"]; Z0 = RF["Z0_ohm"]
    I = P.get("proyecto", {}).get("I_max_uA", 180.0) * 1e-6
    P_haz = P["beam"]["T_MeV"] * 1e6 * I
    P_par = Pcu / RF["factor_Q_real"]
    P_tot = P_par + P_haz
    om = 2 * math.pi * f * 1e6
    rho, dz, Bn = flujo_lazo(g, s, f, P)
    r1, r2 = (x / 1e3 for x in RF["acoplador_rho"])
    sel = (rho >= r1) & (rho <= r2)
    Bn = np.nan_to_num(Bn)
    drho, ddz = rho[1] - rho[0], dz[1] - dz[0]
    phi_cum = np.cumsum(Bn[sel].sum(axis=0) * drho * ddz)               # Φ(h) para lazos de altura h desde la tapa
    fem = om * phi_cum
    Pext = fem ** 2 / (2 * Z0)
    h_crit = float(np.interp(P_tot, Pext, dz)) if Pext[-1] >= P_tot else float("nan")
    Q0 = om * W / P_par; Qext_crit = om * W / P_tot; Qcarg = Q0 * Qext_crit / (Q0 + Qext_crit)
    # sin haz (β > 1 con el lazo crítico a plena corriente): reflexión con el haz apagado
    beta0 = P_tot / P_par
    gamma0 = (beta0 - 1) / (beta0 + 1)
    out = dict(f_MHz=f, Z0_ohm=Z0, P_cobre_ideal_kW=Pcu / 1e3, factor_Q_real=RF["factor_Q_real"],
               P_paredes_kW=P_par / 1e3, I_haz_uA=I * 1e6, P_haz_kW=P_haz / 1e3, P_total_kW=P_tot / 1e3,
               Q0_real=Q0, Q_ext=Qext_crit, Q_cargado=Qcarg, ancho_banda_kHz=f * 1e3 / Qcarg,
               lazo_rho_mm=list(RF["acoplador_rho"]), lazo_h_critico_mm=h_crit * 1e3,
               fem_critica_V=math.sqrt(2 * Z0 * P_tot),
               P_directa_kW=P_tot / 1e3,
               beta_sin_haz=beta0, reflexion_sin_haz=gamma0, P_reflejada_sin_haz_kW=gamma0 ** 2 * P_tot / 1e3,
               Bn_medio_tapa_mT=float(np.mean(Bn[sel][:, :5]) * 1e3))
    log(f"[acoplador] P paredes {P_par / 1e3:.1f} kW + haz {P_haz / 1e3:.2f} kW = {P_tot / 1e3:.1f} kW; "
        f"lazo {RF['acoplador_rho'][0]:.0f}–{RF['acoplador_rho'][1]:.0f} mm → h crítica {h_crit * 1e3:.0f} mm; "
        f"Q0 {Q0:.0f}, Q_carg {Qcarg:.0f}, BW {f * 1e3 / Qcarg:.1f} kHz")
    return out
