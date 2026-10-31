#!/usr/bin/env python3
"""Render del modelo Monte Carlo: búnker cortado, máquina con cuadrante abierto, fuentes de neutrones."""
import sys, pathlib, json
import numpy as np, pyvista as pv, cadquery as cq
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from kiube import load_params, ROOT
from kiube.beam import extraction
from kiube.mc import mc_volumes, target_sources
from build import COLORS
pv.OFF_SCREEN = True
P, _ = load_params(); E = extraction(P)
vols, info = mc_volumes(P, E, "openmc")
S = target_sources(P, E)
COLORS["concrete_portland"] = (0.78, 0.76, 0.72)
Qv = cq.Solid.makeBox(20000, 20000, 20000, pnt=cq.Vector(-10000, -20000, -10000))     # quita y<0 del búnker
Qm = cq.Solid.makeBox(3000, 3000, 4000, pnt=cq.Vector(0, -3000, -2000))              # cuadrante de la máquina
pl = pv.Plotter(off_screen=True, window_size=(1900, 1300)); pl.set_background("white")
for v in vols:
    if v["mat"] in ("helium", "vacuum"):
        continue
    sh = v["shape"].cut(Qv if v["name"] == "vault_concrete" else Qm)
    if sh.Volume() < 1e-6:
        continue
    vv, ff = sh.tessellate(4.0 if v["name"] in ("vault_concrete",) or v["name"].startswith("yoke") else 1.0, 0.25)
    if not ff:
        continue
    m = pv.PolyData(np.array([(p.x, p.y, p.z) for p in vv]), np.hstack([np.full((len(ff), 1), 3), np.array(ff)]).ravel())
    col = (0.25, 0.45, 0.85) if v["name"].startswith("E_") else COLORS.get(v["mat"], (0.7, 0.7, 0.7))
    pl.add_mesh(m, color=col, specular=0.15, opacity=1.0)
for s in S["sources"]:
    pl.add_mesh(pv.Sphere(radius=28, center=s["pos_mm"]), color="#ff2d55")
pl.camera_position = [(5200, -7400, 4300), (0, 300, -200), (0, 0, 1)]
pl.add_text("Modelo Monte Carlo: búnker de concreto (EST), máquina, vacío explícito; rojo = 8 fuentes n (H2-18O)", font_size=13, color="black")
pl.screenshot(ROOT / "out" / "render_mc.png")
print("ok")
