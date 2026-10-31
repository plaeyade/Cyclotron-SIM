"""Vista del CAD de la región central (v0.32): cortes horizontales de las piezas (puente + puntas + mordazas del puller,
chimeneas con su rendija, vástagos) a z = 0, 6 y 10 mm. Entrada: out/print/stl. Salida: validation/cad_region_central.png"""
import pathlib
import numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; D = ROOT / "out/print/stl"
piezas = {"dee_bridge_upper": "#c87533", "dee_bridge_lower": "#c87533", "ion_source_0": "#5b6170", "ion_source_1": "#5b6170"}
M = {n: trimesh.load(D / f"{n}.stl") for n in piezas if (D / f"{n}.stl").exists()}
fig, ax = plt.subplots(1, 3, figsize=(16, 5.6))
for a, z in zip(ax, (0.0, 6.0, 10.0)):
    for n, m in M.items():
        sec = m.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
        if sec is None:
            continue
        p2, _ = sec.to_planar(to_2D=np.eye(4))
        for poly in p2.polygons_full:
            x, y = poly.exterior.xy; a.fill(x, y, color=piezas[n], alpha=0.85, lw=0.4, ec="k")
            for h in poly.interiors:
                xh, yh = h.xy; a.fill(xh, yh, color="w", lw=0.3, ec="k")
    a.set_aspect("equal"); a.set_xlim(-30, 30); a.set_ylim(-30, 30); a.grid(alpha=.3)
    a.set_title(f"z = {z:.0f} mm", fontsize=10); a.set_xlabel("x [mm]")
ax[0].set_ylabel("y [mm]")
fig.suptitle("CAD v0.32 — región central: puntas de los dees + cuello + mordazas del puller con bisel (cobre); punta de la chimenea Ø6 con rendija y vástago (gris)", fontsize=10)
fig.tight_layout(); fig.savefig(ROOT / "validation/cad_region_central.png", dpi=140)
print("ok")
