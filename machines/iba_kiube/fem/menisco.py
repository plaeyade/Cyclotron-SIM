"""Menisco del plasma en la rendija de la fuente PIG interna (v0.30) — cálculo local autoconsistente.

Física
  La fuente H⁻ emite una densidad de corriente j_p fijada por el plasma (Kuo & Hendry, Cyclotrons 1992: fuente
  interna H⁻ con rendija 0.6 × 4 mm, 2.8 mA equivalentes DC, 144 mA/cm²). El límite de Child–Langmuir del gap
  fuente–puller (40 kV, ≈3.2 mm) es ≈4 A/cm²: el plasma está muy por debajo (3–20 % según la fase) y el campo de
  extracción penetra en la rendija. El borde del plasma (menisco) se retira hasta donde la carga espacial de los
  propios iones anula el campo: condición de emisión limitada por carga espacial, ∂φ/∂n = 0 en el menisco.
  Cuasiestático: el tránsito por el gap (≈2 ns) es corto frente al periodo de RF (24.6 ns) → un menisco por nivel
  de voltaje V = V_dee·cos φ.

Método (iteración tipo «frontera de plasma móvil»)
  1. Dominio local alrededor de la rendija (chimenea a tierra, mordazas del puller a V, caja con el potencial del
     modelo global ×V). El plasma detrás del menisco se impone con una penalización de volumen (φ = 0) sobre una
     malla fija: no hay que remallar al mover el menisco.
  2. Se emiten macropartículas desde el menisco con j_p uniforme, se siguen en E + B (1.35 T) y se deposita su carga.
  3. Poisson con esa carga; se mide ∂φ/∂n delante del menisco en nodos de z; donde tira (> 0) el menisco retrocede
     y donde frena avanza. Se ajusta s_m(z) = c0 + c1 z² + c2 z⁴ + c3 z⁶ y se repite hasta |∂φ/∂n| ≪ V/gap.
  Limitaciones: menisco plano a lo ancho de la rendija (0.6 mm; solo se resuelve la forma en z, la que fija la óptica
  vertical); plasma frío (sin temperatura iónica); electrones extraídos no incluidos (el imán los devuelve).
Uso: python fem/menisco.py [--jp 1440] [--cos 0.97 0.86 0.7 0.5]  → out/central_region/menisco/menisco.json
"""
from __future__ import annotations
import argparse, json, math, pathlib, sys, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import central_region as cr

EPS0 = 8.8541878128e-12
QE = 1.602176634e-19
M_H = 1.6735575e-27 + 9.109e-31                 # H⁻
J_P_KUO92 = 1440.0                               # A/m² (144 mA/cm²)
K_CL = 4 * EPS0 / 9 * math.sqrt(2 * QE / M_H)      # Child–Langmuir: j = K φ^1.5 / d²


def frame(info, pr):
    F = info["fuentes"][0]
    S = np.array(F["S"]); en = np.array(F["e_n"]); et = np.array([-en[1], en[0]])
    return S, en, et


def s_m_eval(c, z):
    z2 = np.asarray(z) ** 2
    return c[0] + c[1] * z2 + c[2] * z2 ** 2 + c[3] * z2 ** 3


class Local:
    def __init__(self, pr, info, Phi_glob, h_canal=0.08, h_gap=0.2, liberar=None, rejilla=None):
        from netgen.occ import Box, Cylinder, Pnt, Vec, Z, Axis, OCCGeometry
        from ngsolve import Mesh
        self.pr = pr
        S, en, et = frame(info, pr); self.S, self.en, self.et = S, en, et
        rc, gap, L = pr["chimney_r"], pr["gap_sp"], pr["puller_len"]
        w, h = pr["slit_w"], pr["slit_h"]
        self.s_cav = pr.get("menisco_s_cav", rc - 1.2)          # lo más hondo que puede ir el menisco
        self.s0, self.s1 = rc - 2.2, rc + gap + L + 2.5
        self.uh, self.zh = 4.5, 6.0
        # potencial global (por voltio de dee) en una rejilla que cubre el dominio local → condición en la caja;
        # se calcula ANTES de mallar y luego se libera el modelo global (memoria)
        xs, ys = [], []
        for s in (self.s0, self.s1):
            for u in (-self.uh, self.uh):
                P_ = S + s * en + u * et; xs.append(P_[0]); ys.append(P_[1])
        lo = np.array([min(xs) - 0.5, min(ys) - 0.5, -self.zh - 0.5]); hi = np.array([max(xs) + 0.5, max(ys) + 0.5, self.zh + 0.5])
        n = np.maximum(((hi - lo) / 0.25).astype(int) + 1, 2)
        gx, gy, gz = [np.linspace(lo[i], hi[i], n[i]) for i in range(3)]
        X, Y, Zg = np.meshgrid(gx, gy, gz, indexing="ij")
        if rejilla is not None:
            self._grid = rejilla
        else:
            self._grid = (lo, hi, Phi_glob(np.c_[X.ravel(), Y.ravel(), Zg.ravel()]).reshape(X.shape))
        if liberar:
            liberar()
        import gc; gc.collect()
        ang = math.degrees(math.atan2(en[1], en[0]))
        place = lambda sol: sol.Rotate(Axis((0, 0, 0), Z), ang).Move(Vec(S[0], S[1], 0))
        box = place(Box(Pnt(self.s0, -self.uh, -self.zh), Pnt(self.s1, self.uh, self.zh)))
        for f in box.faces:
            f.name = "caja"
        chim = Cylinder(Pnt(0, 0, -self.zh - 1), Z, r=rc, h=2 * self.zh + 2)
        chim = place(chim) - place(Box(Pnt(self.s_cav, -w / 2, -h / 2), Pnt(rc + 1.0, w / 2, h / 2)))
        for f in chim.faces:
            f.name = "tierra"
        jaws = [place(b.Move(Vec(rc + gap, 0, 0))) for b in cr.mordazas_puller(pr, self.zh + 1)]
        for j in jaws:
            for f in j.faces:
                f.name = "dee"
        air = box - chim - jaws[0] - jaws[1]
        air.mat("vacio")
        from netgen.meshing import MeshingParameters
        mp = MeshingParameters(maxh=0.6, grading=0.3)
        for s in np.arange(self.s_cav - 0.1, rc + 0.6, h_canal * 1.2):
            for u in np.arange(-w / 2 - 0.1, w / 2 + 0.11, h_canal * 1.2):
                for z in np.arange(-h / 2 - 0.3, h / 2 + 0.31, h_canal * 1.2):
                    P_ = S + s * en + u * et
                    mp.RestrictH(float(P_[0]), float(P_[1]), float(z), h_canal)
        for s in np.arange(rc, rc + gap + L + 0.5, 0.3):
            for u in np.arange(-1.5, 1.51, 0.3):
                for z in np.arange(-4.0, 4.01, 0.3):
                    P_ = S + s * en + u * et
                    mp.RestrictH(float(P_[0]), float(P_[1]), float(z), h_gap)
        t0 = time.time()
        self.mesh = Mesh(OCCGeometry(air).GenerateMesh(mp))
        print(f"[menisco] malla local {self.mesh.ne} tetraedros ({time.time() - t0:.0f} s)", flush=True)
        self._setup()

    def _setup(self):
        from ngsolve import H1, L2, GridFunction, x, y, z, Integrate, CoefficientFunction, VoxelCoefficient
        self.fes = H1(self.mesh, order=2, dirichlet="tierra|dee|caja")
        self.gfu = GridFunction(self.fes)
        self.rho = GridFunction(L2(self.mesh, order=0))
        self.vol = np.array(Integrate(CoefficientFunction(1.0), self.mesh, element_wise=True))   # mm³
        S, en, et = self.S, self.en, self.et
        self.s_cf = (x - S[0]) * en[0] + (y - S[1]) * en[1]
        self.u_cf = (x - S[0]) * et[0] + (y - S[1]) * et[1]
        self.z_cf = z
        lo, hi, vals = self._grid
        # VoxelCoefficient espera el arreglo con índices [z][y][x]
        self.phi_g = VoxelCoefficient(tuple(lo), tuple(hi), np.ascontiguousarray(vals.transpose(2, 1, 0)), linear=True)

    def chi(self, c):
        from ngsolve import IfPos
        pr = self.pr
        sm = c[0] + c[1] * self.z_cf ** 2 + c[2] * self.z_cf ** 4 + c[3] * self.z_cf ** 6
        if len(c) > 4:                                               # curvatura a lo ancho de la rendija
            sm = sm + c[4] * (2 * self.u_cf / pr["slit_w"]) ** 2
        inside = IfPos(sm - self.s_cf, 1.0, 0.0) * IfPos(pr["slit_w"] / 2 + 0.05 - self.u_cf, 1.0, 0.0) * \
            IfPos(pr["slit_w"] / 2 + 0.05 + self.u_cf, 1.0, 0.0) * IfPos(pr["slit_h"] / 2 + 0.05 - self.z_cf, 1.0, 0.0) * \
            IfPos(pr["slit_h"] / 2 + 0.05 + self.z_cf, 1.0, 0.0) * IfPos(self.s_cf - self.s_cav + 0.3, 1.0, 0.0)
        return inside

    def solve(self, V, c, con_carga=True):
        from ngsolve import BilinearForm, LinearForm, grad, dx, Preconditioner, TaskManager, solvers
        u, v = self.fes.TnT()
        g = self.gfu
        from ngsolve import BND
        g.vec[:] = 0.0      # un solo Set: cada Set con definedon reinicia el resto del vector
        g.Set(self.mesh.BoundaryCF({"caja": V * self.phi_g, "dee": V, "tierra": 0.0}, default=0.0), BND)
        gam = 2e3                                                    # penalización φ = 0 en el plasma [1/mm²]
        a = BilinearForm(grad(u) * grad(v) * dx + gam * self.chi(c) * u * v * dx)
        f = LinearForm(self.fes)
        if con_carga:
            f += (self.rho / EPS0 * 1e-6) * v * dx                  # ρ [C/m³] / ε0 → V/m² → ×1e-6 V/mm²
        pre = Preconditioner(a, "h1amg")
        with TaskManager():
            a.Assemble(); f.Assemble()
            r = f.vec.CreateVector(); r.data = f.vec - a.mat * g.vec
            inv = solvers.CGSolver(a.mat, pre.mat, maxiter=3000, tol=1e-10, printrates=False)
            du = g.vec.CreateVector(); du.data = inv * r
            fd = np.array(self.fes.FreeDofs(), dtype=bool)
            g.vec.FV().NumPy()[fd] += du.FV().NumPy()[fd]
        return inv.iterations

    def grad_at(self, P):
        from ngsolve import grad
        mp = self.mesh(P[:, 0], P[:, 1], P[:, 2])
        ok = mp["nr"] >= 0
        G = np.zeros((len(P), 3))
        if ok.any():
            G[ok] = np.asarray(grad(self.gfu)(mp[ok])).reshape(-1, 3)
        return G, ok, mp["nr"]

    def phi_at(self, P):
        mp = self.mesh(P[:, 0], P[:, 1], P[:, 2]); ok = mp["nr"] >= 0
        o = np.zeros(len(P))
        if ok.any():
            o[ok] = np.asarray(self.gfu(mp[ok])).ravel()
        return o

    def puntos(self, c, Z_, U_):
        w = self.pr["slit_w"]
        return s_m_eval(c, Z_) + (c[4] if len(c) > 4 else 0.0) * (2 * U_ / w) ** 2

    def j_cl(self, c, Z_, U_, d0=0.1):
        """Child–Langmuir local (emisión limitada por carga espacial): j = K φ(d0)^1.5 / d0², φ a d0 delante del menisco."""
        s_ = self.puntos(c, Z_, U_) + d0
        P = np.c_[self.S[None, :] + s_[:, None] * self.en[None, :] + U_[:, None] * self.et[None, :], Z_]
        ph = np.maximum(self.phi_at(P), 0.0)
        return K_CL * ph ** 1.5 / (d0 * 1e-3) ** 2, ph

    def emitir(self, c, jp, nz=28, nu=5, Bz=1.35, dt_ns=0.002, t_max_ns=10.0, d0=0.1):
        """Emisión limitada por carga espacial desde el menisco: cada haz elemental lleva min(j_p, j_CL) y arranca a d0
        con la energía eφ(d0). Devuelve la carga depositada por elemento [C], el haz de salida y la corriente emitida."""
        pr = self.pr; S, en, et = self.S, self.en, self.et
        w, h = pr["slit_w"], pr["slit_h"]
        zz = (np.arange(nz) + 0.5) / nz * h - h / 2
        uu = (np.arange(nu) + 0.5) / nu * w - w / 2
        Z_, U_ = np.meshgrid(zz, uu, indexing="ij"); Z_ = Z_.ravel(); U_ = U_.ravel()
        jcl, ph = self.j_cl(c, Z_, U_, d0)
        I_k = np.minimum(jp, jcl) * (h / nz * 1e-3) * (w / nu * 1e-3)
        s_ = self.puntos(c, Z_, U_) + d0
        X = np.c_[S[None, :] + s_[:, None] * en[None, :] + U_[:, None] * et[None, :], Z_] * 1e-3
        Vv = np.c_[np.outer(np.sqrt(2 * np.maximum(ph, 1.0) * QE / M_H), en), np.zeros(len(Z_))]
        Q = np.zeros(self.mesh.ne)
        destino = {}
        alive = I_k > 0; dt = dt_ns * 1e-9
        salida = []
        B = np.array([0, 0, Bz])
        for k in range(int(t_max_ns / dt_ns)):
            idx = np.flatnonzero(alive)
            if not len(idx):
                break
            G, ok, nr = self.grad_at(X[idx] * 1e3)
            Q[nr[ok]] += I_k[idx][ok] * dt
            E = -G * 1e3                                                      # V/m
            a_ = (-QE / M_H) * (E + np.cross(Vv[idx], B))
            Vv[idx] += a_ * dt; X[idx] += Vv[idx] * dt
            Xl = X[idx] * 1e3
            s_now = (Xl[:, 0] - S[0]) * en[0] + (Xl[:, 1] - S[1]) * en[1]
            u_now = (Xl[:, 0] - S[0]) * et[0] + (Xl[:, 1] - S[1]) * et[1]
            borde = (s_now > self.s1 - 0.1) | (np.abs(u_now) > self.uh - 0.1) | (np.abs(Xl[:, 2]) > self.zh - 0.1)
            en_rendija = (np.abs(u_now) < w / 2) & (np.abs(Xl[:, 2]) < h / 2) & (s_now < pr["chimney_r"])
            atras = en_rendija & (s_now < self.puntos(c, Xl[:, 2], u_now) - 0.02)   # vuelve al plasma (reflejada)
            fuera = ~ok | borde | atras
            for j in np.flatnonzero(fuera):
                i = idx[j]
                sj = s_now[j]
                dest = ("reflejada" if atras[j] else
                        "salida" if (borde[j] and sj > pr["chimney_r"] + pr["gap_sp"] + 0.5) else
                        "borde_antes_del_gap" if borde[j] else
                        "rendija" if sj < pr["chimney_r"] + 0.05 else
                        "puller" if sj > pr["chimney_r"] + pr["gap_sp"] - 0.05 else "chimenea_exterior")
                destino[dest] = destino.get(dest, 0.0) + float(I_k[i])
                if borde[j] and s_now[j] > pr["chimney_r"] + pr["gap_sp"] + 0.5:      # pasó el gap → haz extraído
                    vs = float(np.dot(Vv[i, :2], en))
                    salida.append(dict(z0=float(Z_[i]), u0=float(U_[i]), z=float(X[i, 2] * 1e3),
                                       zp=float(Vv[i, 2] / max(abs(vs), 1.0)), I=float(I_k[i])))
            alive[idx[fuera]] = False
        vivos = float(I_k[alive].sum())
        if vivos:
            destino["sin_salir_t_max"] = vivos
        return Q, salida, float(I_k.sum()), dict(j_cl=jcl, Z=Z_, U=U_, destino_mA={k: 1e3 * v for k, v in destino.items()})


def ajustar(zk, sk):
    A = np.c_[np.ones_like(zk), zk ** 2, zk ** 4, zk ** 6]
    return np.linalg.lstsq(A, sk, rcond=None)[0]


def menisco_para(L, V, jp, it_max=16, alfa=0.10, log=print):
    """Mueve el menisco hasta que la corriente de Child–Langmuir local iguale a j_p en todo él (s(u, z) = pol(z) + a_u (2u/w)²)."""
    pr = L.pr
    rc, w, h = pr["chimney_r"], pr["slit_w"], pr["slit_h"]
    zk = np.linspace(0.0, 0.92 * h / 2, 10)                                  # nodos en u = 0
    zu = np.linspace(0.0, 0.8 * h / 2, 4); uu = np.full_like(zu, 0.7 * w / 2)  # nodos hacia el borde de la rendija
    c = np.r_[ajustar(zk, np.full_like(zk, pr["plasma_depth"])), 0.0]        # arranque: plano en la cara del plasma
    L.rho.vec[:] = 0.0
    hist = []
    for it in range(it_max):
        L.solve(V, c, con_carga=it > 0)
        Q, sal, I_tot, ex = L.emitir(c, jp)
        rho_new = -Q / (L.vol * 1e-9)                                         # C/m³ (Q = Σ I·dt; H⁻: carga negativa)
        r = L.rho.vec.FV().NumPy()
        r[:] = rho_new if it == 0 else 0.5 * r + 0.5 * rho_new
        L.solve(V, c, con_carga=True)
        j0, _ = L.j_cl(c, zk, np.zeros_like(zk)); ju, _ = L.j_cl(c, zu, uu)
        l0 = np.log(np.maximum(j0, 1e-6 * jp) / jp); lu = np.log(np.maximum(ju, 1e-6 * jp) / jp)
        s0 = L.puntos(c, zk, np.zeros_like(zk)); su = L.puntos(c, zu, uu)
        res = float(max(np.max(np.abs(l0)), np.max(np.abs(lu))))
        Iout = 1e3 * sum(x["I"] for x in sal)
        hist.append(dict(it=it, s_centro=float(s0[0]), s_borde_z=float(s0[-1]), a_u=float(c[4]), max_abs_ln_jcl_jp=res,
                         I_salida_mA=Iout, I_emitida_mA=1e3 * I_tot, I_plasma_mA=1e3 * jp * w * h * 1e-6,
                         destino_mA=ex["destino_mA"]))
        log(f"[menisco] V {V / 1e3:.1f} kV it {it}: s(0) {s0[0]:.3f} s(z borde) {s0[-1]:.3f} a_u {c[4]:.3f} mm  "
            f"max|ln jCL/jp| {res:.2f}  I emitida {1e3 * I_tot:.2f}, salida {Iout:.2f} mA  "
            + " ".join(f"{k} {v:.2f}" for k, v in ex["destino_mA"].items()))
        if it >= 4 and res < 0.15:
            break
        ds0 = np.clip(-alfa * l0, -0.15, 0.15); dsu = np.clip(-alfa * lu, -0.15, 0.15)
        s0n = np.clip(s0 + ds0, L.s_cav + 0.05, rc - 0.03)
        cz = ajustar(zk, s0n)
        sun = np.clip(su + dsu, L.s_cav + 0.05, rc - 0.01)
        au = float(np.clip(np.mean((sun - s_m_eval(cz, zu)) / 0.49), 0.0, 0.6))
        c = np.r_[cz, au]
    zp = [x["zp"] for x in sal]; zz = [x["z"] for x in sal]; Ii = [x["I"] for x in sal]
    return dict(V=V, coef=list(map(float, c)), z_nodos=zk.tolist(), s_nodos=L.puntos(c, zk, np.zeros_like(zk)).tolist(),
                historia=hist, salida=dict(n=len(sal), I_mA=1e3 * sum(Ii),
                                           z_rms_mm=float(np.sqrt(np.average(np.square(zz), weights=Ii))) if sal else None,
                                           zp_rms_mrad=1e3 * float(np.sqrt(np.average(np.square(zp), weights=Ii))) if sal else None,
                                           corr_z_zp=float(np.corrcoef(zz, zp)[0, 1]) if len(sal) > 2 else None))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jp", type=float, default=J_P_KUO92, help="densidad de corriente del plasma [A/m²]")
    ap.add_argument("--cos", type=float, nargs="*", default=[0.97, 0.86, 0.7, 0.5])
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--it", type=int, default=16)
    ap.add_argument("--recalcular", action="store_true", help="recalcula el potencial global (si no, usa la rejilla guardada)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pr = dict(cr.PRESETS["kiube"]); pr["_nombre"] = "kiube"
    for kv in a.set:
        k, v = kv.split("="); pr[k] = type(pr[k])(v)
    out = pathlib.Path(a.out or ROOT / "out/central_region/menisco"); out.mkdir(parents=True, exist_ok=True)
    # potencial global (por voltio de dee) con la cavidad de la rendija profunda y plasma plano → condición de borde local
    pr["menisco_s_cav"] = pr.get("menisco_s_cav", pr["chimney_r"] - 1.2)
    prg = dict(pr, _menisco=[pr["plasma_depth"], 0.0, 0.0, 0.0, 0.0])
    cache = out / "rejilla_phi_global.npz"
    rej = None
    if cache.exists() and "--recalcular" not in sys.argv:
        d = np.load(cache); rej = (d["lo"], d["hi"], d["vals"])
        G = dict(info=cr.build(prg)[1])
    else:
        G = dict(zip(("mesh", "gfu", "info", "st"), cr.solve(prg, None)))
    info = G["info"]
    def Phi(P):
        mesh, gfu = G["mesh"], G["gfu"]
        mp = mesh(P[:, 0], P[:, 1], P[:, 2]); ok = mp["nr"] >= 0
        o = np.zeros(len(P))
        if ok.any():
            o[ok] = np.asarray(gfu(mp[ok])).ravel()
        return o
    L = Local(pr, info, Phi, liberar=lambda: (G.pop("gfu", None), G.pop("mesh", None)), rejilla=rej)
    if rej is None:
        np.savez(cache, lo=L._grid[0], hi=L._grid[1], vals=L._grid[2])
    res = dict(jp_A_m2=a.jp, fuente_jp="Kuo & Hendry, Cyclotrons 1992 (144 mA/cm²)" if a.jp == J_P_KUO92 else "usuario",
               niveles=[])
    for cs in a.cos:
        R = menisco_para(L, pr["V_dee"] * cs, a.jp, it_max=a.it)
        R["cos"] = cs
        res["niveles"].append(R)
        json.dump(res, open(out / f"menisco_jp{int(a.jp)}.json", "w"), indent=1, default=float)
    print(json.dumps([dict(cos=r["cos"], s=r["s_nodos"][::3], salida=r["salida"]) for r in res["niveles"]], default=float))
