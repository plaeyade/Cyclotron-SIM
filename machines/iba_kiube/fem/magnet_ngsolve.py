"""Magnetostática 3D no lineal del imán del KIUBE — misma física que Opera-3d/TOSCA, con NGSolve.

Formulación A (potencial vector, elementos de arista de Nédélec), hierro no lineal con la curva
B-H de materials.yaml (densidad de energía w(B) = ∫H dB, Newton con búsqueda lineal), bobina
homogeneizada con J azimutal, sector de 90° de centro de colina a centro de colina (45°–135°)
con condición PERIÓDICA de rotación y mitad superior (z ≥ 0):
  z = 0          → condición natural n×H = 0  (B normal al plano medio, simetría par)
  aire exterior  → n×A = 0 (flujo tangente)
  φ = 45°/135°   → periódica (rotación de 90° alrededor de z)
Geometría: dominios conformes del LOD fem_mag (out/fem_mag/sector90/domains/*.step), en mm → m.

Uso:  python fem/magnet_ngsolve.py [--maxh_gap 6] [--orden 1] [--NI 105000]
Salidas (out/fem_mag/ngsolve/): mapa Bz(r, φ) en z = 0, ⟨B⟩(r), flutter, comparación con
S09 (B_colina/B_valle medidos) y KEY22 (⟨B⟩ y F del KIUBE), y archivo de mapa para tracking.
"""
from __future__ import annotations

import argparse, json, math, pathlib, sys, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def build_geometry(dom_dir, maxh):
    from netgen.occ import OCCGeometry, Glue, Rotation, Axis, Pnt, Z, IdentificationType
    meta = json.load(open(dom_dir.parent / "domains.json"))
    solids = []
    for d in meta:
        s = OCCGeometry(str(dom_dir / f"{d['name']}.step")).shape
        role = d["role"]
        s.mat("iron" if role == "iron" else "coil" if role == "coil" else "air")
        nm = d["name"]
        if role == "iron":
            fine = nm.startswith(("D_pole", "E_insert", "F_plug"))
            s.maxh = maxh["polo"] if fine else maxh["sector"] if nm.startswith("C_") else maxh["yugo"]
        elif role == "coil":
            s.maxh = maxh["bobina"]
        else:
            s.maxh = maxh["aire"]
        solids.append(s)
    shape = Glue(solids)
    # nombres de caras de borde y refinamiento en la cara de polo
    tol = 0.5
    f45, f135 = [], []
    for f in shape.faces:
        c = f.center
        bb = f.bounding_box
        zs = (bb[0][2], bb[1][2])
        if abs(zs[0]) < tol and abs(zs[1]) < tol:
            f.name = "medio"
            continue
        # plano φ = 45°: x = y ; plano φ = 135°: x = −y (todas las coords de la caja sobre el plano)
        pts = [v.p for v in f.vertices]
        if pts and all(abs(p[0] - p[1]) < tol for p in pts) and c[0] > -tol:
            f.name = "per45"; f45.append(f); continue
        if pts and all(abs(p[0] + p[1]) < tol for p in pts) and c[1] > -tol:
            f.name = "per135"; f135.append(f); continue
    # borde exterior del aire: caras a r = R_aire o en el techo del dominio
    bbx = shape.bounding_box
    zmax = bbx[1][2]
    Rair = max(abs(bbx[0][0]), abs(bbx[1][0]), bbx[1][1])
    for f in shape.faces:
        if f.name in ("medio", "per45", "per135"):
            continue
        pts = np.array([[v.p[0], v.p[1], v.p[2]] for v in f.vertices]) if f.vertices else np.zeros((0, 3))
        fb = f.bounding_box
        if abs(fb[0][2] - zmax) < tol and abs(fb[1][2] - zmax) < tol:
            f.name = "exterior"
        elif len(pts) and np.all(np.hypot(pts[:, 0], pts[:, 1]) > Rair - 1.0):
            f.name = "exterior"
    trafo = Rotation(Axis((0, 0, 0), Z), 90)
    used = set()
    for fa in f45:
        ca = np.array([fa.center[0], fa.center[1], fa.center[2]])
        tgt = np.array([-ca[1], ca[0], ca[2]])                    # rotación +90°
        best, bd = None, 1e9
        for j, fb in enumerate(f135):
            if j in used:
                continue
            cb = np.array([fb.center[0], fb.center[1], fb.center[2]])
            dd = np.linalg.norm(cb - tgt)
            if dd < bd:
                best, bd = j, dd
        if best is None or bd > 1.0:
            raise RuntimeError(f"cara periódica sin pareja (Δ = {bd:.2f} mm)")
        used.add(best)
        fa.Identify(f135[best], "periodica", IdentificationType.PERIODIC, trafo)
    return shape, len(f45)


def bh_energy(materials, curva="bh_acero_bajoC_kiube", Js=None):
    from ngsolve import BSpline
    import yaml
    M = yaml.safe_load(open(ROOT / "materials.yaml"))
    bh = M[curva]
    H = np.array(bh["H"], float)
    B = np.array(bh["B"], float)
    mu0 = 4e-7 * math.pi
    if Js:                                                        # polarización de saturación J_s [T]: B = μ0H + J(H)·J_s/J_s,curva
        J = B - mu0 * H
        B = mu0 * H + J * (Js / J[-1])
    # extrapolación lineal con μ0 más allá del último punto (saturación)
    # v0.33: hasta H ≈ 1e9 A/m (B ≈ 1260 T): un elemento con flujo concentrado en una arista aguda no sale de la tabla
    # (fuera de ella el spline no es convexo: el Newton se estancaba o daba NaN)
    ext = np.array([2e5, 1e6, 5e6, 2e7, 1e8, 1e9])
    Bx = np.concatenate([B, B[-1] + mu0 * ext])
    Hx = np.concatenate([H, H[-1] + ext])
    h_of_b = BSpline(2, [0.0] + list(Bx), list(Hx))
    return h_of_b.Integrate(), (Bx, Hx)


def ampere(mesh, B, bh, P, phi_deg=55.0, n=600):
    """Reparto de ∮H·dl = NI/2 por tramos del camino de flujo (mitad superior, cierra por el plano medio, donde
    H tangencial = 0). Camino a φ = phi_deg (dentro de la colina): r_a en z 0 → z_A (mitad de la tapa) → r_b → z 0."""
    Bb, Hb = bh
    mu0 = 4e-7 * math.pi
    Y, CO = P["yoke"], P["coil"]
    H2 = Y["height"] / 2e3
    zA = H2 - Y["plate_thickness"] / 2e3
    ra = 0.30
    rb = (Y["cavity_radius"] + 0.5 * (Y["side"] / 2 - Y["cavity_radius"])) / 1e3
    c, s_ = math.cos(math.radians(phi_deg)), math.sin(math.radians(phi_deg))
    tramos = [("vertical r=0.30 (gap + polo + sector C + tapa A)", (ra, 0.0), (ra, zA)),
              ("radial en la tapa A", (ra, zA), (rb, zA)),
              ("vertical en el yugo lateral B", (rb, zA), (rb, 0.0))]
    out = []
    mat = mesh.MaterialCF({"iron": 1.0}, default=0.0)
    for nom, (r0, z0), (r1, z1) in tramos:
        t = np.linspace(0, 1, n + 1); tm = 0.5 * (t[1:] + t[:-1])
        dl = np.array([(r1 - r0) * c, (r1 - r0) * s_, z1 - z0]) / n
        hdl_fe = hdl_air = 0.0; Bmax = 0.0
        for tt in tm:
            r = r0 + (r1 - r0) * tt; z = max(z0 + (z1 - z0) * tt, 2e-4)
            mp = mesh(r * c, r * s_, z)
            b = np.array([B[k](mp) for k in range(3)]); bn = float(np.linalg.norm(b))
            if mat(mp) > 0.5:
                h = float(np.interp(bn, Bb, Hb)) if bn > 0 else 0.0
                hdl_fe += (h / bn if bn > 0 else 0.0) * float(b @ dl); Bmax = max(Bmax, bn)
            else:
                hdl_air += float(b @ dl) / mu0
        out.append(dict(tramo=nom, Hdl_hierro=hdl_fe, Hdl_aire=hdl_air, Bmax_hierro=Bmax))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxh_gap", type=float, default=7.0)
    ap.add_argument("--orden", type=int, default=1)
    ap.add_argument("--NI", type=float, nargs="*", default=None,
                    help="uno o varios A·v totales (2 bobinas); por defecto params.yaml. Varios = barrido con arranque en caliente")
    ap.add_argument("--out", default=str(ROOT / "out/fem_mag/ngsolve"))
    ap.add_argument("--curva", default=None, help="curva B-H de materials.yaml (por defecto la del acero del yugo)")
    ap.add_argument("--Js", type=float, default=None, help="reescala la polarización de saturación de la curva [T]")
    ap.add_argument("--diag", action="store_true", help="reparto de Ampère por tramos")
    ap.add_argument("--dom", default=str(ROOT / "out/fem_mag/sector90/domains"), help="carpeta de dominios STEP")
    a = ap.parse_args()
    from ngsolve import (Mesh, HCurl, GridFunction, BilinearForm, CoefficientFunction, curl, dx, sqrt,
                         Variation, x, y, TaskManager, SetNumThreads, Parameter, Periodic, ElementId, VOL)
    from kiube import load_params
    P, _ = load_params()
    CO = P["coil"]
    NIs = a.NI if a.NI else [CO.get("NI_FE_calibrado", CO["NI_total"])]
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    dom = pathlib.Path(a.dom)
    t0 = time.time()
    mh = dict(polo=a.maxh_gap, sector=2.5 * a.maxh_gap, yugo=60.0, bobina=40.0, aire=250.0)
    shape, nper = build_geometry(dom, mh)
    from netgen.occ import OCCGeometry
    ngm = OCCGeometry(shape).GenerateMesh(maxh=250.0, grading=0.35)
    ngm.Scale(1e-3)                                               # mm → m
    mesh = Mesh(ngm)
    print(f"[malla] {mesh.ne} tetraedros, {nper} pares periódicos, {time.time() - t0:.0f} s", flush=True)
    assert "exterior" in set(mesh.GetBoundaries())
    SetNumThreads(2)
    fes = Periodic(HCurl(mesh, order=a.orden, dirichlet="exterior", nograds=True))
    print(f"[fes] {fes.ndof} grados de libertad", flush=True)
    u = fes.TrialFunction()
    mu0 = 4e-7 * math.pi
    if a.curva is None:
        import yaml
        MM = yaml.safe_load(open(ROOT / "materials.yaml"))
        a.curva = next(v["magnetic"] for sec in MM.values() if isinstance(sec, dict)
                       for k, v in sec.items() if k == "steel_1010" and isinstance(v, dict) and "magnetic" in v)
    wB, bhpts = bh_energy(None, a.curva, a.Js)
    w, h = CO["radial_width"] / 1e3, CO["height"] / 1e3
    Jp = Parameter(0.0)
    rr = sqrt(x * x + y * y)
    a_ = BilinearForm(fes, symmetric=True)
    a_ += Variation(wB(sqrt(1e-12 + curl(u) * curl(u))) * dx("iron"))
    a_ += Variation((0.5 / mu0) * curl(u) * curl(u) * dx("air|coil"))
    a_ += Variation(1e-2 * u * u * dx)                            # regularización (gauge)
    a_ += Variation(-CoefficientFunction((-Jp * y / rr, Jp * x / rr, 0)) * u * dx("coil"))
    gfu = GridFunction(fes)
    fd = np.array(fes.FreeDofs(), dtype=bool)
    resumen = []
    def newton(NI, tag=""):
        Jp.Set((NI / 2) / (w * h))                                # A/m² por bobina (mitad superior)
        prev, r0 = 1e300, None
        with TaskManager():
            for it in range(40):
                a_.AssembleLinearization(gfu.vec)
                res = gfu.vec.CreateVector(); a_.Apply(gfu.vec, res)
                res.FV().NumPy()[~fd] = 0.0                       # residuo solo en DOF libres
                rn = res.Norm(); r0 = r0 or rn
                print(f"[newton NI={NI:.0f}{tag}] it {it}  |res| {rn:.3e}  rel {rn / r0:.2e}  ({time.time() - t0:.0f} s)", flush=True)
                if not np.isfinite(rn):
                    return False
                if rn / r0 < 1e-5 or (it > 6 and rn / r0 < 1e-4 and rn > 0.9 * prev):
                    return True
                prev = rn
                inv = a_.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
                du = gfu.vec.CreateVector(); du.data = inv * res
                E0 = a_.Energy(gfu.vec); tau = 1.0
                trial = gfu.vec.CreateVector(); ok = False
                while tau > 1e-4:                                 # búsqueda lineal sobre la energía (NaN = no mejora)
                    trial.data = gfu.vec - tau * du
                    Et = a_.Energy(trial)
                    if np.isfinite(Et) and Et < E0:
                        ok = True; break
                    tau *= 0.5
                del inv
                if not ok:
                    return False
                gfu.vec.data = trial
        return True

    def arranque_lineal(NI, mur=200.0):
        """v0.33: solución lineal con μr efectivo (hierro saturado) como punto de partida del Newton."""
        from ngsolve import BilinearForm as BF, LinearForm as LF, curl as cu
        uu, vv = fes.TnT()
        al = BF(fes, symmetric=True)
        al += (1 / (mu0 * mur)) * cu(uu) * cu(vv) * dx("iron")
        al += (1 / mu0) * cu(uu) * cu(vv) * dx("air|coil")
        al += 2e-2 * uu * vv * dx
        Jv = (NI / 2) / (w * h)
        fl = LF(fes); fl += CoefficientFunction((-Jv * y / rr, Jv * x / rr, 0)) * vv * dx("coil")
        with TaskManager():
            al.Assemble(); fl.Assemble()
            inv = al.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
            gfu.vec.data = inv * fl.vec
            del inv

    for NI in NIs:
        x0 = gfu.vec.CreateVector(); x0.data = gfu.vec
        ok = newton(NI)
        if not ok:
            # v0.33: reintento desde la solución lineal con μr efectivo y, si no, continuación en la corriente
            print("[newton] reintento desde la solución lineal con μr = 200", flush=True)
            arranque_lineal(NI)
            ok = newton(NI, tag=" (arranque lineal)")
        if not ok:
            arranque_lineal(0.5 * NI)
            for fr in (0.5, 0.8, 1.0):
                ok = newton(NI * fr, tag=f" rampa {fr:.1f}")
                if not ok:
                    raise RuntimeError("el Newton no converge ni con arranque lineal ni con continuación en la corriente")
        B = curl(gfu)
        rs = np.arange(0.0, 0.62, 0.005)
        phs = np.arange(45.0, 135.0001, 0.5)
        Bz = np.full((len(rs), len(phs)), np.nan)
        for i, r in enumerate(rs):
            for j, ph in enumerate(phs):
                rx = max(r, 1e-5)
                try:
                    Bz[i, j] = B[2](mesh(rx * math.cos(math.radians(ph)), rx * math.sin(math.radians(ph)), 2e-4))
                except Exception:
                    pass
        tag = f"NI{NI / 1e3:.0f}k"
        np.savez(out / f"bz_midplane_{tag}.npz", r_m=rs, phi_deg=phs, Bz_T=Bz, NI_total=NI)
        # v0.35: mapa extendido para la extracción (hasta la cara del yugo), con control de pertenencia a la malla:
        # fuera de la malla (hierro no mallado o fuera del dominio) → NaN, y el material de cada punto
        rx_ = np.arange(0.40, 0.955, 0.005); px_ = np.arange(45.0, 135.0001, 1.0)
        Bx = np.full((len(rx_), len(px_)), np.nan); mat = np.full((len(rx_), len(px_)), "", dtype=object)
        for i, r in enumerate(rx_):
            for j, ph in enumerate(px_):
                x_, y_ = r * math.cos(math.radians(ph)), r * math.sin(math.radians(ph))
                mp_ = mesh(x_, y_, 2e-4)
                if mp_.nr < 0:
                    continue
                try:
                    Bx[i, j] = B[2](mp_); mat[i, j] = mesh[ElementId(VOL, mp_.nr)].mat
                except Exception:
                    pass
        np.savez(out / f"bz_midplane_ext_{tag}.npz", r_m=rx_, phi_deg=px_, Bz_T=Bx, material=mat.astype(str), NI_total=NI)
        Bm = np.nanmean(Bz, axis=1)
        F = np.nanmean((Bz - Bm[:, None]) ** 2, axis=1) / Bm ** 2
        hill = np.nanmean(Bz[:, np.abs(phs - 90) > 35], axis=1)
        vall = Bz[:, np.argmin(np.abs(phs - 90))]
        r_ = dict(NI_total=NI, ndof=fes.ndof, ne=mesh.ne, maxh_gap_mm=a.maxh_gap, tiempo_s=time.time() - t0,
                  r_m=rs.tolist(), Bavg_T=Bm.tolist(), flutter=F.tolist(), B_colina_T=hill.tolist(), B_valle_T=vall.tolist())
        json.dump(r_, open(out / f"resultado_{tag}.json", "w"), indent=1)
        resumen.append({k: r_[k] for k in ("NI_total", "ndof", "maxh_gap_mm")} |
                       dict(Bavg_045=float(np.interp(0.45, rs, Bm)), F_045=float(np.interp(0.45, rs, F)),
                            Bh_045=float(np.interp(0.45, rs, hill)), Bv_045=float(np.interp(0.45, rs, vall)),
                            Bavg_010=float(np.interp(0.10, rs, Bm))))
        if a.diag:
            dg = ampere(mesh, B, bhpts, P)
            resumen[-1]["ampere"] = dg
            for d_ in dg:
                print(f"[ampere] {d_['tramo']}: hierro {d_['Hdl_hierro']:.0f} A, aire {d_['Hdl_aire']:.0f} A, "
                      f"Bmax {d_['Bmax_hierro']:.2f} T", flush=True)
        resumen[-1]["curva"] = a.curva; resumen[-1]["Js"] = a.Js
        print("[fin]", json.dumps({k: v for k, v in resumen[-1].items() if k != "ampere"}), flush=True)
    json.dump(resumen, open(out / "resumen.json", "w"), indent=1)


if __name__ == "__main__":
    main()
