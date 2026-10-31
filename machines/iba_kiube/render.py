#!/usr/bin/env python3
"""Renders de control: vista isométrica sin mitad superior y vista en planta del plano medio."""
import sys, pathlib, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from kiube import load_params, ROOT
from kiube.beam import extraction
from kiube.geometry import build
from build import COLORS

P, _ = load_params()
E = extraction(P)
parts, info = build(P, E, "neutronics")
out = ROOT / "out"


def _unused(shape, tol=2.0):
    v, f = shape.tessellate(tol, 0.3)
    v = np.array([(p.x, p.y, p.z) for p in v])
    return v[np.array(f)] if f else np.zeros((0, 3, 3))


# ---------- planta del plano medio (sección z=+1 mm y z=+60 mm superpuestas)
import cadquery as cq
orb = E["orbit"]
fig, ax = plt.subplots(figsize=(11, 11), dpi=130)
for zcut, alpha in ((60.0, 0.35), (1.0, 1.0)):
    plane = cq.Face.makePlane(4000, 4000, basePnt=(0, 0, zcut))
    for p in parts:
        if p["name"].endswith("lower"):
            continue
        sec = p["shape"].intersect(plane)
        for f in sec.Faces():
            for w in [f.outerWire()] + f.innerWires():
                pts = np.array([(q.x, q.y) for q in w.discretize(200)]) if hasattr(w, "discretize") else None
                if pts is None:
                    pts = np.array([(q.x, q.y) for q in [w.positionAt(t) for t in np.linspace(0, 1, 200)]])
                col = COLORS.get(p["mat"], (0.7, 0.7, 0.7))
                if p["name"].startswith("E_insert"):
                    col = (0.25, 0.45, 0.85)
                ax.fill(pts[:, 0], pts[:, 1], color=col, alpha=alpha * 0.85, lw=0)
                ax.plot(pts[:, 0], pts[:, 1], color="k", lw=0.3, alpha=alpha)
ax.plot(orb[:, 0], orb[:, 1], color="#e0204a", lw=1.2, label="órbita cerrada 18 MeV (hard-edge)")
for i, pth in enumerate(E["paths"]):
    ax.plot(pth[:, 0], pth[:, 1], color="#18c44a", lw=1.6, label="trayectorias extraídas" if i == 0 else None)
for i, t in enumerate(E["targets"]):
    ax.annotate(f"T{i+1}", (t["x"], t["y"]), xytext=(t["x"] * 1.12, t["y"] * 1.12), fontsize=9, ha="center")
for i, p in enumerate(E["ports"]):
    ax.plot(p["x"], p["y"], "o", ms=7, mfc="none", mec="#d4006a", mew=1.6, label="puerto medido TUP04 Fig.2 (objetivo)" if i == 0 else None)
ax.set_aspect("equal"); ax.set_xlim(-900, 900); ax.set_ylim(-900, 900)
ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
ax.set_title("Sección z = +1 mm (sólido) y z = +60 mm (tenue): polos, insertos, cortes, cámara, blancos")
ax.legend(loc="lower left", fontsize=9)
fig.tight_layout()
fig.savefig(out / "render_midplane.png", bbox_inches="tight", facecolor="white")
print("ok")
