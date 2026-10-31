"""Resonador RF del KIUBE: modos propios del dominio de vacío (fem/rf_dominio.py) con NGSolve (fem/rf_eigen.py).
Identifica el modo en fase (dees con la misma tensión: modo acelerador con h = 2 y dees a 180°) y el de
contrafase; calcula tensión de dee, Q (paredes de cobre) y potencia para 40 kV.
Uso: python fem/rf_kiube.py [--maxh 0.05] [--orden 2] [--param dee.top_z=60 ...]"""
from __future__ import annotations
import argparse, copy, json, math, pathlib, sys, time
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import ngsolve as ng                                                             # noqa: E402
from netgen.occ import OCCGeometry                                               # noqa: E402
from rf_eigen import modos, etiquetar, C0                                        # noqa: E402
from rf_dominio import exportar_step                                             # noqa: E402
from kiube import load_params                                                    # noqa: E402
from kiube.beam import extraction                                                # noqa: E402

MU0, EPS0, SIGMA_CU = 4e-7 * math.pi, 8.8541878128e-12, 5.8e7


def tension_dee(gf, P, j):
    """Tensión aceleradora del dee j: ∫ E·dl en el plano medio (z = 0) a lo largo de un arco r = 0.6·r_max, desde el
    eje del dee (dentro de la apertura, E ≈ 0) hasta el centro de la colina vecina (dentro del gap de colina, E ≈ 0).
    Es la diferencia de potencial dee–tierra que ve el haz (complejo)."""
    DE = P["dee"]; az0 = math.radians(DE["valleys"][j]); r = 0.6 * DE["r_max"] / 1000
    az1 = az0 + math.radians(45.0)                       # centro de la colina siguiente (4 sectores)
    mesh = gf.space.mesh
    n = 120; V = 0
    for i in range(n):
        a = az0 + (i + 0.5) * (az1 - az0) / n
        Ex, Ey, _ = gf(mesh(r * math.cos(a), r * math.sin(a), 1e-4))
        V += (-Ex * math.sin(a) + Ey * math.cos(a)) * r * (az1 - az0) / n
    return V


def resolver(P, E, maxh=0.05, orden=2, f_obj=40.6, log=print, tag="base", n_modos=4, shift_frac=0.6, grading=0.55):
    out = ROOT / "out/fem_rf"; out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    # el dominio (CAD de toda la máquina) se construye en un proceso hijo: la memoria de OpenCascade no vuelve al
    # sistema y la factorización del solver necesita casi todos los 7 GB
    import multiprocessing as mp
    step = out / f"dominio_{tag}.step"
    ctx = mp.get_context("fork")
    pr = ctx.Process(target=exportar_step, args=(P, E, step)); pr.start(); pr.join()
    if pr.exitcode != 0:
        raise RuntimeError(f"exportar_step falló (código {pr.exitcode})")
    sh = OCCGeometry(str(step)).shape
    etiquetar(sh)
    mesh = ng.Mesh(OCCGeometry(sh).GenerateMesh(maxh=maxh, grading=grading))   # 0.55 (v0.19): −33 % de elementos
    mesh.Curve(orden)
    r = modos(sh, n=n_modos, order=orden, mesh=mesh, f_obj=f_obj, shift_frac=shift_frac, log=log)
    res = dict(tag=tag, maxh=maxh, orden=orden, ndof=r["ndof"], tets=mesh.ne, modos=[])
    w = 2 * math.pi
    fs, gs = r["f_MHz"], r["modos"]
    Vs = [[tension_dee(g, P, j) for j in (0, 1)] for g in gs]
    # clasificación: dees en fase (V1/V0 ≈ +1, modo acelerador con h = 2), en contrafase (≈ −1) u otro modo.
    # Si los dos primeros modos están casi degenerados (dees sin acoplar) el solver devuelve mezclas arbitrarias:
    # se construyen las combinaciones con V0 = ±V1.
    cand = []
    if len(fs) >= 2 and abs(fs[1] - fs[0]) < 0.003 * fs[0]:
        a_, b_ = gs[0], gs[1]; Va, Vb = Vs[0], Vs[1]
        wa_tot = 0
        for nombre, sg in (("en fase", 1), ("contrafase", -1)):
            ca, cb = (Vb[0] - sg * Vb[1]), -(Va[0] - sg * Va[1])
            g = ng.GridFunction(a_.space); g.vec.data = ca * a_.vec + cb * b_.vec
            wa, wb = abs(ca) ** 2, abs(cb) ** 2
            cand.append((nombre, (wa * fs[0] + wb * fs[1]) / (wa + wb), g))
        log(f"[rf] {tag}: modos de dee casi degenerados ({fs[0]:.3f} / {fs[1]:.3f} MHz) → combinaciones ±")
    else:
        for f, g, V in zip(fs, gs, Vs):
            if abs(V[0]) < 1e-12 and abs(V[1]) < 1e-12:
                continue
            rr = V[1] / V[0] if abs(V[0]) > 1e-30 else complex(1e9)
            if abs(abs(rr) - 1) < 0.25:
                nombre = "en fase" if rr.real > 0 else "contrafase"
            else:
                nombre = f"otro (|V1/V0| = {abs(rr):.2f})"
            if any(c[0] == nombre for c in cand):
                continue
            cand.append((nombre, f, g))
    for nombre, f, g in cand:
        V = [tension_dee(g, P, j) for j in (0, 1)]
        s = 40e3 / max(abs(V[0]), 1e-30)
        om = w * f * 1e6
        We = 2 * 0.25 * EPS0 * ng.Integrate(ng.InnerProduct(g, g).real, mesh) * s ** 2         # energía eléctrica media
        Hc = ng.BoundaryFromVolumeCF(ng.curl(g)) / (om * MU0)
        Rs = math.sqrt(om * MU0 / (2 * SIGMA_CU))
        Ploss = 2 * 0.5 * Rs * ng.Integrate(ng.InnerProduct(Hc, Hc).real * ng.ds("pec"), mesh) * s ** 2
        Q = om * 2 * We / Ploss if Ploss > 0 else float("nan")
        res["modos"].append(dict(f_MHz=f, fase=nombre, V1_sobre_V0=abs(V[1]) / max(abs(V[0]), 1e-30),
                                 Q_cobre=Q, P_kW_40kV=Ploss / 1e3, W_J_40kV=2 * We))
        if nombre == "en fase":
            res["campo_en_fase"] = g; res["escala_40kV"] = s
    res["modos"].sort(key=lambda m: (m["fase"] != "en fase", m["f_MHz"]))
    res["f_bruto_MHz"] = fs
    res["t_s"] = round(time.time() - t0, 1)
    log(f"[rf] {tag}: " + "; ".join(f"{m['f_MHz']:.3f} MHz ({m['fase']}, Q {m['Q_cobre']:.0f}, {m['P_kW_40kV']:.1f} kW)" for m in res["modos"]))
    return res, r


def figura(g, P, ruta, res, n=140):
    """|E| (normalizado al máximo) en el plano medio y en el corte vertical por el eje del dee 0."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mesh = g.space.mesh

    def ev(x, y, z):
        try:
            mp = mesh(x, y, z)
            if mp.nr < 0:
                return np.nan
            v = g(mp)
            return math.sqrt(sum(abs(c) ** 2 for c in v))
        except Exception:                                        # noqa: BLE001
            return np.nan
    R = P["vacuum_chamber"]["inner_radius"] / 1000
    xs = np.linspace(-R, R, n)
    A = np.array([[ev(x, y, 5e-4) for x in xs] for y in xs])
    az = math.radians(P["dee"]["valleys"][0])
    H = P["yoke"]["height"] / 2000; zt = H + (P["exterior"]["rf_housing_h"] - 20) / 1000
    rs = np.linspace(0, R, n); zs = np.linspace(0, zt, n)
    B = np.array([[ev(r * math.cos(az), r * math.sin(az), z) for r in rs] for z in zs])
    m = np.nanmax(np.concatenate([A.ravel(), B.ravel()]))
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
    im = ax[0].imshow(A / m, origin="lower", extent=[-R, R, -R, R], cmap="inferno", vmin=0, vmax=1)
    ax[0].set_title("|E| en el plano medio (z = 0)"); ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("y [m]")
    ax[1].imshow(B / m, origin="lower", extent=[0, R, 0, zt], cmap="inferno", vmin=0, vmax=1, aspect="auto")
    ax[1].set_title("|E| en el corte por el eje del dee (stem a la derecha)"); ax[1].set_xlabel("r [m]"); ax[1].set_ylabel("z [m]")
    fig.colorbar(im, ax=ax, shrink=0.8, label="|E| / máx")
    mf = res["modos"][0]
    fig.suptitle(f"KIUBE — modo acelerador (dees en fase): {mf['f_MHz']:.3f} MHz, Q(Cu) = {mf['Q_cobre']:.0f}, "
                 f"{mf['P_kW_40kV']:.1f} kW a 40 kV")
    fig.savefig(ruta, dpi=140)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxh", type=float, default=0.05); ap.add_argument("--orden", type=int, default=2)
    ap.add_argument("--param", nargs="*", default=[]); ap.add_argument("--tag", default="base")
    ap.add_argument("--figura", default=None, help="PNG con |E| del modo en fase (plano medio y corte por el dee)")
    a = ap.parse_args()
    P, _ = load_params(); P = copy.deepcopy(P)
    for kv in a.param:
        k, v = kv.split("="); g, n = k.split(".")
        P[g][n] = float(v)
    E = extraction(P)
    res, _ = resolver(P, E, a.maxh, a.orden, tag=a.tag, n_modos=3, shift_frac=0.5)
    if "rf" in P and "campo_en_fase" in res:
        from rf_acoplos import acoplador
        res["acoplador"] = acoplador(res, res["campo_en_fase"], res["escala_40kV"], P)
    if a.figura:
        figura(res["campo_en_fase"], P, a.figura, res)
    res.pop("campo_en_fase", None); res.pop("escala_40kV", None)
    json.dump(res, open(ROOT / f"out/fem_rf/resultado_{a.tag}.json", "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))
