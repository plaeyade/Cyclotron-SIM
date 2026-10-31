"""Modos propios de RF del resonador del Cyclone KIUBE (NGSolve, elementos de Nédélec).

Problema: rot rot E = k² E en el vacío entre conductores (PEC: n × E = 0). Plano medio z = 0: pared magnética
(condición natural), válida para el modo acelerador (E tangencial en el plano medio, simétrico arriba/abajo).
Autovalores con Arnoldi (desplazamiento-inversión) cerca de la frecuencia buscada; el núcleo de gradientes (λ ≈ 0)
se descarta. f = c·√λ / 2π (geometría en metros).

Uso como módulo:  modos(shape_m, n=4, maxh=..., order=2) → frecuencias [MHz] y autovectores
Verificación:     python fem/rf_eigen.py --verificar   (cavidad cilíndrica TM010 y coaxial λ/4)
"""
from __future__ import annotations

import math, sys, time
from netgen.occ import OCCGeometry, Cylinder, Box, Pnt, Axes, X, Y, Z  # noqa: F401
import ngsolve as ng

C0 = 299792458.0


def etiquetar(shape, z_sym=0.0, tol=1e-6):
    """Caras en el plano medio → 'sim' (pared magnética, natural); resto → 'pec'."""
    n = 0
    for f in shape.faces:
        (x0, y0, z0), (x1, y1, z1) = f.bounding_box
        if abs(z0 - z_sym) < tol and abs(z1 - z_sym) < tol:
            f.name = "sim"; n += 1
        else:
            f.name = "pec"
    return shape


def modos(shape, n=4, maxh=0.05, order=2, minh=None, grading=0.3, log=print, mesh=None, f_obj=40.0, shift_frac=0.8, inverse="sparsecholesky"):
    """shape en metros, con caras ya etiquetadas ('pec' / 'sim'). Busca los n modos más cercanos a f_obj [MHz]
    (Arnoldi con desplazamiento-inversión, campo complejo). Devuelve dict con f_MHz, malla y autovectores."""
    t0 = time.time()
    if mesh is None:
        geo = OCCGeometry(shape)
        mp = dict(maxh=maxh, grading=grading)
        if minh:
            mp["minh"] = minh
        mesh = ng.Mesh(geo.GenerateMesh(**mp))
        mesh.Curve(order)
    log(f"[rf] malla: {mesh.ne} tetraedros, contornos {sorted(set(mesh.GetBoundaries()))} ({time.time() - t0:.0f} s)")
    fes = ng.HCurl(mesh, order=order, dirichlet="pec", complex=True)
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.curl(u) * ng.curl(v) * ng.dx)
    m = ng.BilinearForm(u * v * ng.dx)
    gf = ng.GridFunction(fes, multidim=n + 1)
    k0 = 2 * math.pi * f_obj * 1e6 / C0
    with ng.TaskManager():
        a.Assemble(); m.Assemble()
        # Arnoldi con desplazamiento-inversión cerca de la frecuencia buscada; el núcleo de gradientes (λ = 0)
        # aparece como un autovalor ~0 y se descarta
        lam = ng.ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), list(gf.vecs), shift=k0 * k0 * shift_frac, inverse=inverse)
    pares = sorted(((C0 * math.sqrt(abs(l.real)) / (2 * math.pi) / 1e6, i) for i, l in enumerate(lam)))
    pares = [(fr, i) for fr, i in pares if fr > 0.05 * f_obj][:n]
    f = [fr for fr, _ in pares]
    log(f"[rf] {fes.ndof} gdl, f = {', '.join(f'{x:.3f}' for x in f)} MHz ({time.time() - t0:.0f} s)")
    vecs = [gf.vecs[i] for _, i in pares]
    gfs = []
    for vv in vecs:
        g = ng.GridFunction(fes); g.vec.data = vv; gfs.append(g)
    return dict(f_MHz=f, mesh=mesh, fes=fes, modos=gfs, ndof=fes.ndof)


def verificar(log=print):
    """Casos con solución exacta (mitad z ≥ 0, pared magnética en z = 0)."""
    out = {}
    # 1) cavidad cilíndrica cerrada (todo PEC) R = 0.3 m, h = 0.1 m: TM010 = 2.405 c / (2π R)
    R, h = 0.3, 0.1
    cyl = Cylinder(Pnt(0, 0, -h / 2), Z, r=R, h=h)
    for fc in cyl.faces:
        fc.name = "pec"
    r = modos(cyl, n=2, maxh=0.06, order=3, log=log, f_obj=380.0)
    ex = 2.404826 * C0 / (2 * math.pi * R) / 1e6
    out["TM010"] = dict(fem=r["f_MHz"][0], exacto=ex, error_rel=r["f_MHz"][0] / ex - 1)
    # Q con paredes de cobre (σ = 5.8e7 S/m): Q = (η/Rs) · 2.405 / (2 (1 + R/h))
    gf, mesh = r["modos"][0], r["mesh"]
    om = 2 * math.pi * r["f_MHz"][0] * 1e6; mu0 = 4e-7 * math.pi; Rs = math.sqrt(om * mu0 / (2 * 5.8e7))
    Hc = ng.BoundaryFromVolumeCF(ng.curl(gf)) / (om * mu0)
    P = 0.5 * Rs * ng.Integrate(ng.InnerProduct(Hc, Hc).real * ng.ds("pec"), mesh)
    W = 0.5 * 8.8541878128e-12 * ng.Integrate(ng.InnerProduct(gf, gf).real, mesh)
    Qex = 376.730313 / Rs * 2.404826 / (2 * (1 + R / h))
    out["Q_TM010_cobre"] = dict(fem=om * W / P, exacto=Qex, error_rel=om * W / P / Qex - 1)
    # 2) coaxial λ/4: conductor interior cortocircuitado arriba (z = L) y extremo abierto en z = 0 (pared magnética)
    #    modo TEM: f = c / (4 L)
    L, a_, b_ = 0.75, 0.02, 0.06
    outer = Cylinder(Pnt(0, 0, 0), Z, r=b_, h=L)
    inner = Cylinder(Pnt(0, 0, 0), Z, r=a_, h=L)
    coax = etiquetar(outer - inner)
    r = modos(coax, n=2, maxh=0.05, order=3, log=log, f_obj=100.0)
    ex = C0 / (4 * L) / 1e6
    out["coax_lambda4"] = dict(fem=r["f_MHz"][0], exacto=ex, error_rel=r["f_MHz"][0] / ex - 1)
    return out


if __name__ == "__main__":
    if "--verificar" in sys.argv:
        import json
        res = verificar()
        import pathlib
        out = pathlib.Path(__file__).resolve().parents[1] / "validation/rf_verificacion.json"
        json.dump(res, open(out, "w"), indent=1)
        print(json.dumps(res, indent=1))
