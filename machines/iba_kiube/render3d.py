#!/usr/bin/env python3
"""Renders 3D (pyvista offscreen): máquina abierta y corte en cuadrante."""
import sys, pathlib
import numpy as np
import pyvista as pv

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from kiube import load_params, ROOT
from kiube.beam import extraction
from kiube.geometry import build
from build import COLORS

pv.OFF_SCREEN = True
P, _ = load_params()
E = extraction(P)
parts, info = build(P, E, "neutronics")
out = ROOT / "out"


def mesh(shape, tol):
    v, f = shape.tessellate(tol, 0.2)
    if not f:
        return None
    v = np.array([(p.x, p.y, p.z) for p in v])
    f = np.hstack([np.full((len(f), 1), 3), np.array(f)]).ravel()
    return pv.PolyData(v, f)


meshes = []
for p in parts:
    if p["mat"] == "helium":
        continue
    tol = 3.0 if p["name"].startswith(("yoke", "leg", "lifting")) else 0.8
    m = mesh(p["shape"], tol)
    if m is None:
        continue
    col = COLORS.get(p["mat"], (0.7, 0.7, 0.7))
    if p["name"].startswith("E_insert"):
        col = (0.25, 0.45, 0.85)
    meshes.append((p["name"], m, col))


def lines(pl):
    o = E["orbit"]
    pl.add_mesh(pv.lines_from_points(np.c_[o, np.full(len(o), 0.5)]), color="#e0204a", line_width=3)
    for pth in E["paths"]:
        pl.add_mesh(pv.lines_from_points(np.c_[pth, np.full(len(pth), 0.5)]), color="#18c44a", line_width=4)


# 1) máquina abierta: todo salvo mitad superior
pl = pv.Plotter(off_screen=True, window_size=(1800, 1400))
pl.set_background("white")
for n, m, c in meshes:
    if n.endswith("upper") or n.startswith(("leg", "lifting")):
        continue
    pl.add_mesh(m, color=c, smooth_shading=False, specular=0.2)
lines(pl)
pl.camera_position = [(1900, -2200, 1900), (0, 0, -100), (0, 0, 1)]
pl.add_text("Cyclone KIUBE — mitad inferior abierta (LOD neutrónico)", font_size=14, color="black")
pl.screenshot(out / "render_open.png")

# 2) máquina completa con corte de un cuadrante
pl = pv.Plotter(off_screen=True, window_size=(1800, 1400))
pl.set_background("white")
import cadquery as cq
Q = cq.Solid.makeBox(3000, 3000, 4000, pnt=cq.Vector(0, -3000, -2000))
for p in parts:
    if p["mat"] == "helium":
        continue
    sh = p["shape"].cut(Q)
    if sh.Volume() < 1e-6:
        continue
    tol = 3.0 if p["name"].startswith(("yoke", "leg", "lifting")) else 0.8
    m = mesh(sh, tol)
    if m is None:
        continue
    col = COLORS.get(p["mat"], (0.7, 0.7, 0.7))
    if p["name"].startswith("E_insert"):
        col = (0.25, 0.45, 0.85)
    pl.add_mesh(m, color=col, specular=0.2)
lines(pl)
pl.camera_position = [(2100, -2400, 1100), (0, 0, -100), (0, 0, 1)]
pl.add_text("Corte de cuadrante — yugo, sectores, polos, bobinas, cámara, dees, blancos", font_size=14, color="black")
pl.screenshot(out / "render_cut.png")

# 3) modelo FEM magnético: sector 90° con aire translúcido
from kiube.geometry import build_fem_mag
doms = build_fem_mag(P, E, sector=True)
pl = pv.Plotter(off_screen=True, window_size=(1800, 1400))
pl.set_background("white")
for p in doms:
    m = mesh(p["shape"], 6.0 if p["name"] in ("air", "yoke_AB_upper") else 1.0)
    if m is None:
        continue
    if p["name"] == "air":
        pl.add_mesh(m.extract_feature_edges(25), color=(0.25, 0.45, 0.85), line_width=2)
    else:
        col = (0.25, 0.45, 0.85) if p["name"].startswith("E_") else COLORS.get(p["mat"], (0.7, 0.7, 0.7))
        pl.add_mesh(m, color=col, specular=0.2)
pl.camera_position = [(3800, -2600, 3200), (0, 1100, 500), (0, 0, 1)]
pl.add_text("fem_mag · sector 90° (colina 45°–colina 135°), mitad superior; azul = contorno del dominio de aire", font_size=14, color="black")
pl.screenshot(out / "render_fem_mag.png")
print("ok")
