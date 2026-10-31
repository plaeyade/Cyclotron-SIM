"""Ajuste fino de isocronía (v0.34): jacobiano medido sobre las órbitas de equilibrio.

Perillas: Δα por lado [°] en R2 = 100…510 mm (10 radios, interpolación lineal) y el contorno del escalón 2 L [mm].
Base: el mapa de la iteración 3 con L = 520 (misma geometría con esta parametrización).
  python fem/shim_iso2.py jac i        → caso con la perilla i perturbada (i = 1…10: Δα; 11: L)
  python fem/shim_iso2.py paso TAG_BASE → resuelve δ con el jacobiano y calcula el caso nuevo (TAG_BASE → TAG_BASE+)
Objetivo: Δf/f(E) = c (nivel libre) en 1–18 MeV; peso 2 en 14–18 MeV; Tikhonov sobre δ.
"""
from __future__ import annotations
import json, pathlib, subprocess, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import shim_iso as S1

R2 = [60, 100, 150, 200, 280, 360, 400, 430, 455, 480, 500]      # el nodo de 60 mm queda fijo en 0
VAR = S1.VAR
PASO = [None] + [0.7] * 10 + [4.0]                             # índice 1…10: Δα (0.7° por lado); 11: L (+4 mm)
TS = S1.TS


def cotas_dee(holgura_mm=5.0, L=(515.0, 526.0)):
    """Δα máximo por lado para que el borde de la colina quede a ≥ holgura del borde de la caja del dee (39° en los valles):
    45° + (22.5° − asin(10/r)) + Δα ≤ 70.5° − holgura/r."""
    import math
    mx = [(70.5 - math.degrees(holgura_mm / r)) - 45 - (22.5 - math.degrees(math.asin(10 / r))) for r in R2[1:]]
    return dict(da_min=[-3.0] * 10, da_max=mx, L=L)


def de_viejo(d_old, r_old=S1.R_CTRL):
    return [float(np.interp(r, r_old, d_old)) for r in R2]


BASE_DA = de_viejo([0.0, 2.217, 0.722, 0.074, -0.242, 0.13, -0.621, 3.642])
BASE_L = 520.0


def caso(tag, da, L):
    d = VAR / tag
    npz = d / "ngsolve/bz_midplane_NI115k.npz"
    if not npz.exists():
        sets = [f"pole.shim_r={json.dumps(R2)}", f"pole.shim_dalpha={json.dumps([round(float(x), 4) for x in da])}",
                f"pole.step2_edge_L={L}"]
        subprocess.run([sys.executable, str(ROOT / "fem/fem_mag_variante.py"), "--tag", tag, "--set", *sets], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
        for mh in ("7.0", "6.5", "7.5", "6.0"):
            with open(d / "fe.log", "w") as lg:
                rc = subprocess.run([sys.executable, str(ROOT / "fem/magnet_ngsolve.py"), "--dom", str(d / "sector90/domains"),
                                     "--out", str(d / "ngsolve"), "--maxh_gap", mh], cwd=ROOT, stdout=lg,
                                    stderr=subprocess.STDOUT).returncode
            if rc == 0 and npz.exists():
                break
        else:
            raise RuntimeError(f"FE sin converger en {tag}")
    f = d / "orbitas.json"
    if not f.exists():
        E = S1.orbitas(npz)
        json.dump(dict(tag=tag, r_ctrl=R2, dalpha=list(map(float, da)), L=L, orbitas=E), open(f, "w"), indent=1)
    return json.load(open(f))


def errores(D, Ts=None):
    m = {o["T_MeV"]: o["error_iso_ppm"] for o in D["orbitas"]}
    return m


def jacobiano():
    b = errores(json.load(open(VAR / "shim_it3_L520/orbitas.json")))
    cols = []
    for i in range(1, 12):
        f = VAR / f"j2_{i}" / "orbitas.json"
        if not f.exists():
            return None, None
        e = errores(json.load(open(f)))
        cols.append({T: (e[T] - b[T]) / PASO[i] for T in b if T in e})
    return b, cols


def resolver(b, cols, base_da, base_L, reg=600.0, Ts=None, peso_fase=1.0, regL=1.0, regc=1.0, cotas=None):
    """Minimiza el error punto a punto Δf/f(E) y la FASE ACUMULADA del ion central φ(E) = 720°·Σ Δf/f por vuelta
    (vueltas ≈ (E − E0)/ΔE, ΔE = 4·40 kV·sin 39° ≈ 0.1007 MeV), con nivel libre c, Tikhonov y suavidad del perfil."""
    Ts = Ts or [T for T in sorted(b) if 0.5 <= T <= 18.0 and all(T in c for c in cols)]
    T = np.array(Ts)
    B = np.array([b[t] for t in Ts]) * 1e-6; J = np.array([[c[t] for c in cols] for t in Ts]) * 1e-6
    n = J.shape[1]
    # integración de la fase en una rejilla fina de energía (interpolación lineal del error entre energías)
    Tf = np.arange(0.75, min(18.0, T[-1]) + 1e-9, 0.1007)
    Iw = np.array([np.interp(Tf, T, np.eye(len(T))[k]) for k in range(len(T))]).T     # error en Tf = Iw @ err(T)
    C = np.tril(np.ones((len(Tf), len(Tf)))) * 720.0                                  # fase [°] = C @ err(Tf)
    Cs = C @ Iw
    ones = np.ones(len(T))
    # incógnitas x = [δ (n), c]; err(T) = B + J δ + c
    nf = len(Tf)
    # incógnitas x = [δ (n), c, φ_ini]; fase(Tf) = φ_ini + Cs·err
    A1 = np.c_[J, ones, np.zeros(len(T))] * 1e4                          # 1 = 100 ppm
    A2 = np.c_[Cs @ J, Cs @ ones, np.ones(nf)] / 20.0 * peso_fase        # 1 = 20° de fase
    b1 = -B * 1e4; b2 = -(Cs @ B) / 20.0 * peso_fase
    # Tikhonov: δα (1 = 1°), δL (regL por mm: el escalón 2 no puede pasar del escalón 1, 529 mm) y nivel c
    # (regc por 0.1 %: más campo exige más corriente con el hierro saturado)
    Rg = np.c_[np.diag([reg / 600.0] * (n - 1) + [regL]), np.zeros((n, 2))]
    Rc = np.zeros((1, n + 2)); Rc[0, n] = regc * 1e3
    D2 = np.zeros((n - 3, n + 2))
    for k in range(n - 3):
        D2[k, k:k + 3] = [1, -2, 1]
    A = np.vstack([A1, A2, Rg, 2 * reg / 600.0 * D2, Rc]); bb = np.r_[b1, b2, np.zeros(n), np.zeros(n - 3), 0.0]
    if cotas is None:
        x, *_ = np.linalg.lstsq(A, bb, rcond=None)
    else:                                                               # cotas físicas sobre Δα absoluto y L absoluto
        from scipy.optimize import lsq_linear
        lo = np.full(n + 2, -np.inf); hi = np.full(n + 2, np.inf)
        for i in range(10):
            lo[i] = cotas["da_min"][i] - base_da[i + 1]; hi[i] = cotas["da_max"][i] - base_da[i + 1]
        lo[10] = cotas["L"][0] - base_L; hi[10] = cotas["L"][1] - base_L
        x = lsq_linear(A, bb, bounds=(lo, hi)).x
    d, c, p0 = x[:n], x[n], x[n + 1]
    err = B + J @ d + c
    fase = p0 + Cs @ err
    return np.r_[0.0, np.array(base_da[1:]) + d[:10]], base_L + d[10], c, dict(
        T=Ts, pred_ppm=(err * 1e6).tolist(), fase_pred_max_deg=float(np.abs(fase).max()), delta=d.tolist(), nivel_ppm=float(c * 1e6), fase_ini_deg=float(p0))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "jac":
        i = int(sys.argv[2])
        da = list(BASE_DA); L = BASE_L
        if 1 <= i <= 10:
            da[i] += PASO[i]
        else:
            L += PASO[i]
        D = caso(f"j2_{i}", da, L)
        print(f"[j2_{i}]", [(o["T_MeV"], round(o["error_iso_ppm"])) for o in D["orbitas"]][-6:], flush=True)
    elif cmd == "resolver":
        b, cols = jacobiano()
        kw = dict(reg=float(sys.argv[2]) if len(sys.argv) > 2 else 600.0)
        if len(sys.argv) > 3:
            kw.update(regL=float(sys.argv[3]), regc=float(sys.argv[4]))
        da, L, c, info = resolver(b, cols, BASE_DA, BASE_L, **kw)
        print(json.dumps(dict(dalpha=[round(x, 3) for x in da], L=round(L, 2), **info), default=float))
