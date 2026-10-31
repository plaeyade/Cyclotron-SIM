"""Ajuste de isocronía del imán (v0.33): perfil Δα(r) de los bordes de colina (shimming) y nivel de campo.

Pasos (cada caso FE ≈ 10 min, sector de 90°, no lineal, 115 kA·v):
  python fem/shim_iso.py caso TAG '[Δα_1, …, Δα_8]'   → FE + órbitas de equilibrio de TAG
  python fem/shim_iso.py respuesta                   → base + 1 caso por radio de control (+1° por lado)
  python fem/shim_iso.py ajustar [iteración]          → mínimos cuadrados sobre la respuesta lineal de ⟨B⟩(r)
Objetivo: f_rev(E) = f_RF/h en 1–17.5 MeV, con un factor de escala libre (corriente de la bobina).
Resultados: out/fem_mag_variantes/<TAG>/{ngsolve/bz_midplane_NI115k.npz, orbitas.json}
"""
from __future__ import annotations
import json, math, pathlib, subprocess, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
R_CTRL = [60, 100, 150, 200, 280, 360, 430, 500]                # mm
VAR = ROOT / "out/fem_mag_variantes"
TS = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 17.5, 18]
F_RF, H = 40.640, 2


def fe(tag, dalpha, extra=()):
    d = VAR / tag
    npz = d / "ngsolve/bz_midplane_NI115k.npz"
    if not npz.exists():
        sets = [f"pole.shim_r={json.dumps(R_CTRL)}", f"pole.shim_dalpha={json.dumps([float(x) for x in dalpha])}", *extra]
        subprocess.run([sys.executable, str(ROOT / "fem/fem_mag_variante.py"), "--tag", tag, "--set", *sets], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)
        for mh in ("7.0", "6.5", "7.5"):                            # otra malla si un elemento degenerado da NaN
            with open(d / "fe.log", "w") as lg:
                rc = subprocess.run([sys.executable, str(ROOT / "fem/magnet_ngsolve.py"), "--dom", str(d / "sector90/domains"),
                                     "--out", str(d / "ngsolve"), "--maxh_gap", mh], cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT).returncode
            if rc == 0 and npz.exists():
                break
        else:
            raise RuntimeError(f"FE sin converger en {tag}")
    return npz


def orbitas(npz, escala=1.0, Ts=TS):
    import orbitas as O
    M = O.Mapa(npz, escala=escala)
    out = []; rg = 0.05
    for T in Ts:
        try:
            o = O.orbita_cerrada(M, T, rg)
            if not np.isfinite(o["r0_m"]):
                break
            rg = o["r0_m"]; o["error_iso_ppm"] = 1e6 * (o["f_rev_MHz"] / (F_RF / H) - 1); out.append(o)
        except Exception:
            break
    return out


def bavg(npz):
    d = np.load(npz); return d["r_m"], d["Bz_T"].mean(1)


def caso(tag, dalpha, extra=()):
    npz = fe(tag, dalpha, extra)
    E = orbitas(npz)
    json.dump(dict(tag=tag, dalpha=list(map(float, dalpha)), r_ctrl=R_CTRL, extra=list(extra), orbitas=E),
              open(VAR / tag / "orbitas.json", "w"), indent=1)
    for o in E:
        print(f"[{tag}] T {o['T_MeV']:5.1f}  r {o['r_medio_m'] * 1e3:6.1f}  Δf/f {o['error_iso_ppm']:+7.0f} ppm  νr {o['nu_r']:.3f}  νz {o['nu_z']:.3f}", flush=True)
    return E


def sensibilidad():
    """S_i [T/°] de cada perilla, medida en el FE (respuesta en su radio efectivo, media ±10 mm). La perilla de 60 mm se
    descarta (el Newton no converge: arista fina junto al tapón); la de 500 mm actúa en 470 mm (más fuera, borde)."""
    rb, Bb = bavg(VAR / "shim_base" / "ngsolve/bz_midplane_NI115k.npz")
    r_ef, S = [], []
    for i, rc in enumerate(R_CTRL):
        f = VAR / f"shim_d{i}" / "ngsolve/bz_midplane_NI115k.npz"
        if i == 0 or not f.exists():
            continue
        r, B = bavg(f); d = B - Bb
        re = 0.47 if rc == 500 else rc / 1e3
        m = (r > re - 0.0101) & (r < re + 0.0101)
        r_ef.append(re); S.append(float(d[m].mean()))
    return np.array(r_ef), np.array(S)


def ajustar(base_tag, dalpha0, reg=0.05, Tmax=17.5, fijas=()):
    """Δα_nuevo(r_ctrl) = Δα0 + δ, con ΔB(r) = S(r)·δ(r) ≈ −(Δf/f)(r)·⟨B⟩(r) + c·⟨B⟩(r) en las órbitas de 1…Tmax MeV."""
    r_ef, S = sensibilidad()
    rr_ctrl = np.array(R_CTRL[1:], float) / 1e3                       # perillas activas (100…500 mm)
    r0, B0 = bavg(VAR / base_tag / "ngsolve/bz_midplane_NI115k.npz")
    E = json.load(open(VAR / base_tag / "orbitas.json"))["orbitas"]
    sel = [o for o in E if 1.0 <= o["T_MeV"] <= Tmax]
    ro = np.array([o["r_medio_m"] for o in sel]); err = np.array([o["error_iso_ppm"] for o in sel]) * 1e-6
    # el nivel medio lo corrige la corriente de la bobina: se ajusta solo la forma (mediana de 3–14 MeV restada)
    c0 = float(np.median([o["error_iso_ppm"] for o in sel if 3 <= o["T_MeV"] <= 14])) * 1e-6
    err = err - c0
    Bo = np.interp(ro, r0, B0)
    So = np.interp(ro, r_ef, S)
    n = len(rr_ctrl)
    # δ(r) = Σ_i δ_i hat_i(r) (interpolación lineal entre perillas, constante fuera)
    H_ = np.array([np.interp(ro, rr_ctrl, np.eye(n)[i]) for i in range(n)]).T
    A = np.c_[So[:, None] * H_, -Bo]                                 # incógnitas: δ_1..δ_n, c
    # el nivel c se fija en 0 (columna anulada): nivel y Δα uniforme son casi degenerados; la corriente de la bobina se
    # reajusta después con las órbitas (escala del mapa)
    A[:, -1] = 0.0
    b = -err * Bo
    D = np.zeros((n - 2, n + 1))
    for k in range(n - 2):
        D[k, k:k + 3] = [1, -2, 1]
    w = 1.0 / np.maximum(np.abs(Bo) * 1e-3, 1e-6)                     # pesos: errores relativos de 0.1 % → 1
    Areg = np.vstack([A * w[:, None], reg * 30 * D, reg * np.c_[np.eye(n), np.zeros((n, 1))]])
    breg = np.r_[b * w, np.zeros(n - 2), np.zeros(n)]
    for i in fijas:                                                  # perillas fijas (índice en R_CTRL): δ_i = 0
        Areg[:, i - 1] = 0.0; Areg[len(b) + (n - 2) + i - 1, i - 1] = 1e3
    x, *_ = np.linalg.lstsq(Areg, breg, rcond=None)
    delta, c = x[:n], x[n]
    nuevo = np.r_[0.0, np.array(dalpha0[1:]) + delta]
    pred = err + (A[:, :n] @ delta) / Bo - c
    return nuevo, c, dict(r_orb_mm=(ro * 1e3).tolist(), err_ppm=(err * 1e6).tolist(), pred_ppm=(pred * 1e6).tolist(),
                          delta=delta.tolist(), nivel_c=float(c), S_T_por_grado=S.tolist(), r_ef_m=r_ef.tolist())


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "caso":
        caso(sys.argv[2], json.loads(sys.argv[3]), sys.argv[4:])
    elif cmd == "respuesta":
        i = int(sys.argv[2])
        if i < 0:
            caso("shim_base", [0.0] * len(R_CTRL))
        else:
            d = [0.0] * len(R_CTRL); d[i] = 1.0
            caso(f"shim_d{i}", d)
    elif cmd == "ajustar":
        base = sys.argv[2]; d0 = json.loads(sys.argv[3])
        fijas = tuple(int(x) for x in sys.argv[4].split(",")) if len(sys.argv) > 4 else ()
        nuevo, c, info = ajustar(base, d0, Tmax=12.0 if fijas else 17.5, fijas=fijas)
        print(json.dumps(dict(nuevo=[round(x, 3) for x in nuevo.tolist()], nivel=c, **info), default=float))
