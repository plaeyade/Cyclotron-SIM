"""Órbitas de equilibrio y aceleración completa en el mapa magnético FE del KIUBE (v0.33).

Mapa: Bz(r, φ) en el plano medio (out/fem_mag/ngsolve/bz_midplane_*.npz, 0–0.615 m, periodo 90°). Se descompone en
armónicos 4k (k = 0…K) con splines en r (derivadas analíticas, sin ruido de malla). Fuera del plano medio, expansión de
primer orden: B_r = z ∂Bz/∂r, B_θ = (z/r) ∂Bz/∂θ.

Ecuaciones con θ como variable independiente (Gordon), momentos normalizados al momento total p (p_θ = √(1−p_r²−p_z²)):
  r' = r p_r/p_θ          p_r' = p_θ − (r/Bρ)(B_z − (p_z/p_θ) B_θ)
  z' = r p_z/p_θ          p_z' = (r/Bρ)(B_r − (p_r/p_θ) B_θ)
  t' = r/(v p_θ)
Órbita cerrada: (r, p_r) periódicos en un periodo de 90° (Newton). Sintonías: matriz de transferencia lineal de un
periodo (diferencias finitas) → ν = 4·acos(traza/2)/2π. Frecuencia de revolución: 4·∫dt sobre el periodo.
Uso: python fem/orbitas.py → out/orbitas/equilibrio.json
"""
from __future__ import annotations
import json, math, pathlib, sys
import numpy as np
from scipy.interpolate import make_smoothing_spline
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
C = 299792458.0
M_H = 939.294           # MeV (H⁻)


def brho(T, m=M_H):
    return math.sqrt(T * T + 2 * T * m) / 299.792458        # T·m


class Mapa:
    def __init__(self, f=None, K=10, escala=1.0, lam=None, datos=None):
        if datos is None:
            f = f or ROOT / "out/fem_mag/ngsolve/bz_midplane_NI115k.npz"
            datos = np.load(f)
        d = datos
        self.r = d["r_m"]; ph = np.radians(d["phi_deg"]); B = d["Bz_T"]
        self.escala = escala
        # armónicos 4k sobre el periodo [45°, 135°) (último punto repetido)
        th = ph[:-1]; Bp = B[:, :-1]
        n = Bp.shape[1]
        F = np.fft.rfft(Bp, axis=1) / n
        self.K = K
        # la rejilla cubre un periodo de 90° → el armónico m de la FFT es 4m en θ
        self.a = [2 * F[:, m].real if m else F[:, 0].real for m in range(K + 1)]
        self.b = [-2 * F[:, m].imag if m else 0 * F[:, 0].real for m in range(K + 1)]
        self.th0 = th[0]
        r = self.r
        self.sa = [make_smoothing_spline(r, c, lam=lam) for c in self.a]
        self.sb = [make_smoothing_spline(r, c, lam=lam) for c in self.b]

    def B(self, r, th):
        """Bz, ∂Bz/∂r, ∂Bz/∂θ (T, T/m, T/rad) en arrays."""
        r = np.asarray(r, float); th = np.asarray(th, float)
        Bz = np.zeros_like(r); dr = np.zeros_like(r); dt = np.zeros_like(r)
        for m in range(self.K + 1):
            w = 4 * m; x = w * (th - self.th0)
            c, s = np.cos(x), np.sin(x)
            A, Bb = self.sa[m](r), self.sb[m](r)
            Ad, Bd = self.sa[m](r, 1), self.sb[m](r, 1)
            Bz += A * c + Bb * s
            dr += Ad * c + Bd * s
            dt += w * (-A * s + Bb * c)
        k = self.escala
        return k * Bz, k * dr, k * dt


def mapa_key22(f_fe=None):
    """Mapa sintético con el ⟨B⟩(r) y el flutter F(r) MEDIDOS del KIUBE (KEY22 Fig. 3) y la forma azimutal del mapa FE
    (bordes de colina reales) reescalada al flutter medido: Bz = ⟨B⟩_K [1 + (B_FE/⟨B_FE⟩ − 1)·√(F_K/F_FE)]."""
    d = np.load(f_fe or ROOT / "out/fem_mag/ngsolve/bz_midplane_NI115k.npz")
    r, ph, B = d["r_m"], d["phi_deg"], d["Bz_T"]
    k = np.genfromtxt(ROOT / "validation/kiube_Bavg_flutter_KEY22.csv", delimiter=",", skip_header=4)
    rb = k[~np.isnan(k[:, 1])]; rf = k[~np.isnan(k[:, 2])]
    Bk = np.interp(r, rb[:, 0], rb[:, 1] / 1e4)
    Bk[r < rb[0, 0]] = rb[0, 1] / 1e4
    Fk = np.interp(r, rf[:, 0], rf[:, 2])
    Bm = B.mean(1, keepdims=True); mod = B / Bm - 1
    Ffe = (mod ** 2).mean(1, keepdims=True)
    esc = np.sqrt(np.clip(Fk[:, None], 0, None) / np.maximum(Ffe, 1e-6))
    Bs = Bk[:, None] * (1 + mod * esc)
    sel = r <= 0.50                                            # KEY22 ⟨B⟩ medido hasta 0.497 m
    return dict(r_m=r[sel], phi_deg=ph, Bz_T=Bs[sel])


def rhs(th, Y, M, br):
    r, pr, z, pz = Y[0], Y[1], Y[2], Y[3]
    pt = np.sqrt(np.maximum(1 - pr * pr - pz * pz, 1e-12))
    Bz, dBr, dBt = M.B(r, th)
    Br = z * dBr; Bt = z / r * dBt
    return np.array([r * pr / pt,
                     pt - (r / br) * (Bz - pz / pt * Bt),
                     r * pz / pt,
                     (r / br) * (Br - pr / pt * Bt)]), pt


def rk4_periodo(M, br, Y0, th0, n=360, dth=None, v=None):
    """Integra un periodo de 90° (o dth) con RK4; devuelve Y final, tiempo [s] y la trayectoria."""
    dth = dth or math.pi / 2
    h = dth / n; Y = np.array(Y0, float); th = th0; t = 0.0; tr = [Y.copy()]
    for _ in range(n):
        k1, p1 = rhs(th, Y, M, br); k2, p2 = rhs(th + h / 2, Y + h / 2 * k1, M, br)
        k3, p3 = rhs(th + h / 2, Y + h / 2 * k2, M, br); k4, p4 = rhs(th + h, Y + h * k3, M, br)
        if v:
            t += h / 6 * (Y[0] / p1 + 2 * (Y[0] + h / 2 * k1[0]) / p2 + 2 * (Y[0] + h / 2 * k2[0]) / p3 + (Y[0] + h * k3[0]) / p4) / v
        Y = Y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); th += h
        tr.append(Y.copy())
    return Y, t, np.array(tr)


def orbita_cerrada(M, T, r_guess, th0=None):
    br = brho(T); th0 = M.th0 if th0 is None else th0
    x = np.array([r_guess, 0.0])
    for _ in range(30):
        Yf, _, _ = rk4_periodo(M, br, [x[0], x[1], 0, 0], th0)
        F = Yf[:2] - x
        if np.max(np.abs(F)) < 1e-10:
            break
        J = np.zeros((2, 2)); e = [1e-6, 1e-7]
        for j in range(2):
            xp = x.copy(); xp[j] += e[j]
            Yp, _, _ = rk4_periodo(M, br, [xp[0], xp[1], 0, 0], th0)
            J[:, j] = ((Yp[:2] - xp) - F) / e[j]
        x = x - np.linalg.solve(J, F)
    g = 1 + T / M_H; v = C * math.sqrt(1 - 1 / g ** 2)
    Yf, tq, tr = rk4_periodo(M, br, [x[0], x[1], 0, 0], th0, v=v)
    # matrices de transferencia lineales (r, p_r) y (z, p_z)
    Mr = np.zeros((2, 2)); Mz = np.zeros((2, 2)); e = 1e-6
    for j in range(2):
        d = np.zeros(4); d[j] = e
        Yp, _, _ = rk4_periodo(M, br, np.r_[x, 0, 0] + d, th0); Mr[:, j] = (Yp[:2] - Yf[:2]) / e
        d = np.zeros(4); d[2 + j] = e
        Yp, _, _ = rk4_periodo(M, br, np.r_[x, 0, 0] + d, th0); Mz[:, j] = (Yp[2:] - Yf[2:]) / e
    nu = lambda Mx: 4 * math.acos(max(-1, min(1, np.trace(Mx) / 2))) / (2 * math.pi)
    r_media = float(np.mean(tr[:, 0]))
    return dict(T_MeV=T, r0_m=float(x[0]), pr0=float(x[1]), r_medio_m=r_media, f_rev_MHz=1 / (4 * tq) / 1e6,
                nu_r=nu(Mr), nu_z=nu(Mz), traza_r=float(np.trace(Mr)), traza_z=float(np.trace(Mz)),
                estable_z=bool(abs(np.trace(Mz)) < 2), estable_r=bool(abs(np.trace(Mr)) < 2))


def equilibrio(M, Ts, f_rf_MHz=40.640, h=2):
    out = []; rg = 0.05
    for T in Ts:
        o = orbita_cerrada(M, T, rg); rg = o["r0_m"]
        o["error_iso_ppm"] = 1e6 * (o["f_rev_MHz"] / (f_rf_MHz / h) - 1)
        out.append(o)
        print(f"[orb] T {T:6.2f} MeV  r {o['r_medio_m'] * 1e3:6.1f} mm  f_rev {o['f_rev_MHz']:.4f} MHz  "
              f"Δf/f {o['error_iso_ppm']:+8.0f} ppm  νr {o['nu_r']:.3f}  νz {o['nu_z']:.3f}", flush=True)
    return out


if __name__ == "__main__":
    M = Mapa()
    Ts = [0.5, 1, 2, 3, 4, 6, 8, 10, 12, 14, 15, 16, 17, 18, 18.5]
    E = equilibrio(M, Ts)
    d = ROOT / "out/orbitas"; d.mkdir(parents=True, exist_ok=True)
    json.dump(dict(mapa="out/fem_mag/ngsolve/bz_midplane_NI115k.npz", escala=M.escala, f_rf_MHz=40.640, h=2, orbitas=E),
              open(d / "equilibrio.json", "w"), indent=1)
