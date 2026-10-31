"""Diagnóstico de las pérdidas verticales de la región central (v0.29): dónde (r, azimut), cuándo (vuelta) y con qué
fase y z inicial se pierden las partículas por |z|, y cómo evoluciona z(t) de las que sobreviven.
Uso: python fem/central_region_diag_z.py [clave=valor ...]  → out/central_region/diag_z/{diag.json, diag_z.png}
"""
import json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import central_region as cr
from kiube import load_params

if __name__ == "__main__":
    pr = dict(cr.PRESETS["kiube"]); pr["_nombre"] = "kiube"
    for kv in sys.argv[1:]:
        k, v = kv.split("="); pr[k] = type(pr[k])(v)
    pr["puller_slit_h"] = float(pr.get("puller_slit_h", 10.0))
    out = ROOT / "out/central_region/diag_z"; out.mkdir(parents=True, exist_ok=True)
    P, _ = load_params()
    mesh, gfu, info, st = cr.solve(pr, None)
    Ef = cr.FEField(mesh, gfu, False)
    F = info["fuentes"][0]
    r, trs, B, Bf = cr.run_species(Ef, pr, F, P, (-70.0, 30.0, 5.0), 6.0, 400, 0.9 * pr["R_dom"])
    T_rf = 1e9 / pr["f_rf"]; T_rev = 2 * T_rf                          # h = 2
    filas = []
    for b in B:
        tr = b["tr"]
        if not len(tr):
            continue
        t, x, y, z, T = tr[-1]
        filas.append(dict(motivo=b["motivo"], fase=b["fase"], z0=float(tr[0, 3]), r_fin=float(math.hypot(x, y)),
                          az_fin=float(math.degrees(math.atan2(y, x)) % 360), vuelta=float(t / T_rev), T_keV=float(T),
                          zmax=float(np.abs(tr[:, 3]).max())))
    zl = [f for f in filas if f["motivo"] == "z"]; vv = [f for f in filas if f["motivo"] == "vive"]
    res = dict(n=len(filas), perdidas_z=len(zl), vivas=len(vv),
               vuelta_perdida_z=np.percentile([f["vuelta"] for f in zl], [10, 50, 90]).tolist() if zl else None,
               r_perdida_z_mm=np.percentile([f["r_fin"] for f in zl], [10, 50, 90]).tolist() if zl else None,
               fase_z=np.percentile([f["fase"] for f in zl], [10, 50, 90]).tolist() if zl else None,
               fase_vivas=np.percentile([f["fase"] for f in vv], [10, 50, 90]).tolist() if vv else None,
               z0_abs_z=float(np.mean([abs(f["z0"]) for f in zl])) if zl else None,
               z0_abs_vivas=float(np.mean([abs(f["z0"]) for f in vv])) if vv else None, filas=filas)
    json.dump(res, open(out / "diag.json", "w"), indent=1, default=float)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for f, c in ((zl, "#d0182a"), (vv, "#16a085")):
        ax[0].scatter([g["fase"] for g in f], [abs(g["z0"]) for g in f], s=8, color=c)
    ax[0].set_xlabel("fase de salida [°RF]"); ax[0].set_ylabel("|z| inicial [mm]"); ax[0].set_title("rojo: perdida por |z|; verde: vive")
    ax[1].hist([g["vuelta"] for g in zl], bins=24, color="#d0182a"); ax[1].set_xlabel("vuelta de la pérdida"); ax[1].set_title("¿Cuándo?")
    ax[2].scatter([g["r_fin"] * math.cos(math.radians(g["az_fin"])) for g in zl], [g["r_fin"] * math.sin(math.radians(g["az_fin"])) for g in zl], s=6, color="#d0182a")
    ax[2].set_aspect("equal"); ax[2].set_title("¿Dónde? (plano medio, mm)"); ax[2].grid(alpha=.3)
    for b in B[:60]:
        if b["motivo"] == "vive" and len(b["tr"]):
            pass
    fig.tight_layout(); fig.savefig(out / "diag_z.png", dpi=130)
    # z(t) de algunas vivas y perdidas
    fig, ax = plt.subplots(figsize=(8, 3.6))
    nz = nv = 0
    for b in B:
        if b["motivo"] == "z" and nz < 15:
            ax.plot(b["tr"][:, 0] / T_rev, b["tr"][:, 3], color="#d0182a", lw=0.6); nz += 1
        elif b["motivo"] == "vive" and nv < 15:
            ax.plot(b["tr"][:, 0] / T_rev, b["tr"][:, 3], color="#16a085", lw=0.6); nv += 1
    ax.set_xlabel("vueltas"); ax.set_ylabel("z [mm]"); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(out / "z_t.png", dpi=130)
    print({k: v for k, v in res.items() if k != "filas"})
