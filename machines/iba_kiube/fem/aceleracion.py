"""Aceleración completa del haz H⁻ hasta 18 MeV en el mapa magnético FE (v0.33).

- Campo: fem/orbitas.Mapa (armónicos 4k + splines; B_r = z ∂Bz/∂r, B_θ = (z/r) ∂Bz/∂θ), con un factor de escala (corriente
  de la bobina) elegido para minimizar la deriva de fase.
- RF: 2 dees en fase (h = 2) con cajas de 39° centradas en 90° y 270°; golpe de energía en cada borde: al entrar en el dee
  ΔT = +e·V(t), al salir ΔT = −e·V(t) (H⁻), V(t) = V_dee cos(ωt + φ0). Sin focalización eléctrica de los gaps (despreciable
  por encima de ~1 MeV) ni campo de RF fuera de los bordes.
- Fase: φ = fase RF en el cruce del centro del dee − 90° (φ = 0 → máxima ganancia: 4 eV sin(h·39°/2) por vuelta).
- Arranque: final de la región central (T0 ≈ 0.75 MeV), sobre la órbita de equilibrio con amplitudes radial y vertical
  tomadas de fem/central_region (descentrado ≈ 6 mm; |z| ≤ 2 mm) y fases repartidas en la ventana aceptada.
- Pérdidas: |z| > 11.8 mm (gap de colina ±12), frenado (T < 0.8 T0) o > 400 vueltas. Llegada: T ≥ 18 MeV.
Uso: python fem/aceleracion.py [--npz ...] [--escala auto|x] [--n 200]  → out/orbitas/aceleracion_<tag>.json
"""
from __future__ import annotations
import argparse, json, math, pathlib, sys, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import orbitas as O

C = 299792458.0; QE = 1.602176634e-19; MEV = 1e6 * QE
M_KG = O.M_H * MEV / C ** 2
DEES = (90.0, 270.0); ANCHO = 39.0


def gamma_v(T):
    g = 1 + T / O.M_H
    return g, C * np.sqrt(1 - 1 / g ** 2)


def B3(M, X):
    """B (T) en puntos X [m] (N×3) con la expansión de primer orden en z."""
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    r = np.hypot(x, y) + 1e-12; th = np.arctan2(y, x)
    Bz, dBr, dBt = M.B(r, th)
    Br = z * dBr; Bt = z / r * dBt
    c, s = x / r, y / r
    return np.c_[Br * c - Bt * s, Br * s + Bt * c, Bz]


def dentro_dee(th_deg):
    d = np.zeros_like(th_deg, bool)
    for a in DEES:
        d |= np.abs((th_deg - a + 180) % 360 - 180) < ANCHO / 2
    return d


def acelerar(M, T0, fases, dr0, z0, V_kV=40.0, f_rf=40.640e6, h=2, T_fin=18.0, dt_ns=0.04, max_vueltas=400, semilla=1, log=None,
             estado=None, lamina=None, r_max=0.60):
    """estado (v0.35): dict(X [m], P [kg m/s], T [MeV], phi0 [rad]) por partícula, p. ej. el final de la región central;
    sustituye el arranque sobre la órbita de equilibrio (T0, dr0, z0 se ignoran; fases solo se registra).
    lamina (v0.35): dict(phi_deg, r_in, r_out [m]) = lámina extractora; al cruzarla la partícula se detiene con motivo
    "lamina" y su estado queda en el resultado (estados_lamina); T_fin deja de aplicarse. r > r_max → motivo "r"."""
    rng = np.random.default_rng(semilla)
    N = len(fases)
    w = 2 * math.pi * f_rf; dt = dt_ns * 1e-9; q = -QE
    if estado is None:
        o = O.orbita_cerrada(M, T0, 0.09)
        r0, pr0 = o["r0_m"], o["pr0"]
        th0 = M.th0                                                     # 45° (centro de colina)
        # betatrón: desplazamiento radial (centrado) y vertical aleatorios
        ar = rng.uniform(-1, 1, N) * dr0; az = rng.uniform(-1, 1, N) * z0
        r = r0 + ar
        g, v = gamma_v(T0); p = g * M_KG * v
        pr = pr0 * p; pt = math.sqrt(max(p * p - pr * pr, 0))
        X = np.c_[r * math.cos(th0), r * math.sin(th0), az]
        er = np.array([math.cos(th0), math.sin(th0), 0]); et = np.array([-math.sin(th0), math.cos(th0), 0])
        P = np.tile(pr * er + pt * et, (N, 1))                          # momento [kg m/s] (sentido antihorario)
        # φ0: el RF en el instante t se escribe cos(w t + φ0); se elige φ0 por partícula para que su fase relativa inicial
        # (fase RF en el centro del próximo dee − 90°) sea la pedida → se calibra con el tiempo de llegada al primer centro
        phi0 = np.zeros(N); T = np.full(N, T0); cal = np.zeros(N, bool)
    else:
        X = np.array(estado["X"], float); P = np.array(estado["P"], float); T = np.array(estado["T"], float)
        phi0 = np.array(estado["phi0"], float); cal = np.ones(N, bool); r0 = float(np.hypot(X[:, 0], X[:, 1]).mean())
        T0 = float(T.mean())
    vivo = np.ones(N, bool); llego = np.zeros(N, bool); motivo = np.array(["vive"] * N, dtype=object)
    th_prev = np.degrees(np.arctan2(X[:, 1], X[:, 0])); in_prev = dentro_dee(th_prev)
    vuelta = np.zeros(N); hist = [[] for _ in range(N)]          # (vueltas, T, r, fase, z)
    est_lam = {}
    t = 0.0; nmax = int(max_vueltas * h / f_rf / dt)
    t0 = time.time()
    for k in range(nmax):
        idx = np.flatnonzero(vivo)
        if not len(idx):
            break
        Xi = X[idx]; Pi = P[idx]
        B = B3(M, Xi)
        gm = np.sqrt(1 + np.einsum("ij,ij->i", Pi, Pi) / (M_KG * C) ** 2)
        tv = (q * dt / (2 * M_KG * gm))[:, None] * B
        pp = Pi + np.cross(Pi, tv)
        sv = 2 * tv / (1 + np.einsum("ij,ij->i", tv, tv))[:, None]
        Pi = Pi + np.cross(pp, sv)
        gm = np.sqrt(1 + np.einsum("ij,ij->i", Pi, Pi) / (M_KG * C) ** 2)
        Xi = Xi + Pi / (M_KG * gm)[:, None] * dt
        t += dt
        th = np.degrees(np.arctan2(Xi[:, 1], Xi[:, 0])); ins = dentro_dee(th)
        # cruce del centro de un dee: calibración de φ0 (primera vez) y registro de fase
        for a in DEES:
            dprev = (th_prev[idx] - a + 180) % 360 - 180; dnow = (th - a + 180) % 360 - 180
            cruz = (dprev < 0) & (dnow >= 0) & (np.abs(dnow) < 10)
            for j in np.flatnonzero(cruz):
                i = idx[j]
                if not cal[i]:
                    phi0[i] = math.radians(90.0 + fases[i]) - w * t; cal[i] = True
                fase = (math.degrees(w * t + phi0[i]) - 90.0 + 180) % 360 - 180
                vuelta[i] += 0.5
                rr = float(np.hypot(Xi[j, 0], Xi[j, 1]))
                hist[i].append((vuelta[i], float(T[i]), rr, fase, float(Xi[j, 2])))
        # golpes de RF en los bordes (solo tras calibrar la fase)
        ent = ins & ~in_prev[idx]; sal = ~ins & in_prev[idx]
        for sel, sgn in ((ent, +1.0), (sal, -1.0)):
            for j in np.flatnonzero(sel & cal[idx]):
                i = idx[j]
                dT = sgn * V_kV * 1e-3 * math.cos(w * t + phi0[i])
                Tn = T[i] + dT
                gn, vn = gamma_v(Tn); pn = gn * M_KG * vn
                Pi[j] *= pn / np.linalg.norm(Pi[j]); T[i] = Tn
        th_old = th_prev[idx].copy()
        X[idx] = Xi; P[idx] = Pi; th_prev[idx] = th; in_prev[idx] = ins
        # pérdidas y llegada
        zl = np.abs(Xi[:, 2]) > 0.0118
        fr = T[idx] < 0.8 * T0
        rr_ = np.hypot(Xi[:, 0], Xi[:, 1])
        rm = rr_ > r_max
        if lamina is None:
            fin = T[idx] >= T_fin
        else:                                                       # cruce de la lámina (azimut fijo, tramo radial)
            a = lamina["phi_deg"]
            dprev = (th_old - a + 180) % 360 - 180; dnow = (th - a + 180) % 360 - 180
            fin = (dprev < 0) & (dnow >= 0) & (np.abs(dnow) < 10) & (rr_ >= lamina["r_in"]) & (rr_ <= lamina["r_out"])
            for j in np.flatnonzero(fin):
                est_lam[int(idx[j])] = dict(X=Xi[j].tolist(), P=Pi[j].tolist(), T=float(T[idx[j]]), t=t,
                                            fase=float((math.degrees(w * t + phi0[idx[j]]) - 90.0 + 180) % 360 - 180))
        for j in np.flatnonzero(zl | fr | fin | rm):
            i = idx[j]; vivo[i] = False
            motivo[i] = ("llega" if lamina is None else "lamina") if fin[j] else ("z" if zl[j] else "r" if rm[j] else "frenada")
            llego[i] = fin[j]
        if log and k % 20000 == 0:
            log(f"[acel] t {t * 1e9:7.0f} ns  vivas {vivo.sum()}  llegan {llego.sum()}  T medio {T[vivo].mean() if vivo.any() else 0:.2f} MeV  ({time.time() - t0:.0f} s)")
    motivo[vivo] = "max_vueltas"
    return dict(n=N, llegan=int(llego.sum()), motivos={m: int((motivo == m).sum()) for m in set(motivo)},
                fases=list(map(float, fases)), motivo=list(motivo), vueltas=vuelta.tolist(),
                hist=[h_[::4] + h_[-1:] for h_ in hist], r0_m=r0, T0=T0, estados_lamina=est_lam)


def deriva_fase(E, V_kV=40.0, h=2, phi_ini=0.0, T0=0.75, T1=18.0):
    """Fase del ion central (sin betatrón) por integración vuelta a vuelta: Δφ = 360 h (f_RF/h / f_rev − 1),
    ΔT = 4 V sin(h·39°/2) cos φ. E: órbitas de equilibrio (lista de dicts con T_MeV, f_rev_MHz)."""
    Ts = np.array([o["T_MeV"] for o in E]); fr = np.array([o["f_rev_MHz"] for o in E])
    dE = 4 * V_kV * 1e-3 * math.sin(math.radians(h * ANCHO / 2))
    T, ph, out = T0, phi_ini, []
    for n in range(2000):
        f = np.interp(T, Ts, fr)
        out.append((n, T, ph))
        T += dE * math.cos(math.radians(ph))
        ph += 360 * h * (40.640 / h / f - 1)
        if T >= T1 or T <= 0.5 * T0 or abs(ph) > 90:
            out.append((n + 1, T, ph)); break
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=str(ROOT / "out/fem_mag/ngsolve/bz_midplane_NI115k.npz"))
    ap.add_argument("--escala", default="1.0")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--tag", default="base")
    ap.add_argument("--T0", type=float, default=0.75)
    ap.add_argument("--fases", type=float, nargs=2, default=[-40.0, 40.0])
    a = ap.parse_args()
    M = O.Mapa(a.npz, escala=float(a.escala))
    fases = np.random.default_rng(3).uniform(a.fases[0], a.fases[1], a.n)
    R = acelerar(M, a.T0, fases, dr0=0.006, z0=0.002, log=print)
    R.update(npz=a.npz, escala=float(a.escala))
    d = ROOT / "out/orbitas"; d.mkdir(parents=True, exist_ok=True)
    json.dump(R, open(d / f"aceleracion_{a.tag}.json", "w"), default=float)
    print(json.dumps({k: v for k, v in R.items() if k in ("n", "llegan", "motivos")}))


def mejor_nivel(E, T1=18.0, phis=np.arange(-60, 61, 5), ks_ppm=np.arange(-8000, 8001, 100)):
    """Factor de nivel k (f_rev → k·f_rev, ≈ escala del campo) y fase inicial que llevan al ion central más lejos y,
    entre los que llegan a T1, con la menor excursión de fase. Devuelve (k, φ0, T alcanzada, |φ|máx, vueltas)."""
    mejor = None
    for kp in ks_ppm:
        E2 = [dict(o, f_rev_MHz=o["f_rev_MHz"] * (1 + kp * 1e-6)) for o in E]
        for p0 in phis:
            h_ = deriva_fase(E2, phi_ini=float(p0), T1=T1)
            Tf = h_[-1][1]; ex = max(abs(x[2]) for x in h_)
            clave = (Tf >= T1 - 1e-6, -ex if Tf >= T1 - 1e-6 else Tf)
            if mejor is None or clave > mejor[0]:
                mejor = (clave, 1 + kp * 1e-6, float(p0), Tf, ex, len(h_))
    return mejor[1:]
