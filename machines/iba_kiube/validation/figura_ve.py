"""Figura de la energía variable: planta del puerto VE (#8) con el haz de 18 y de 13 MeV,
la horquilla larga y la brida en cuña. Uso: python validation/figura_ve.py [T_MeV]"""
import math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pyvista as pv, cadquery as cq                                   # noqa: E402
from kiube import load_params                                          # noqa: E402
from kiube.beam import extraction                                      # noqa: E402
from kiube.geometry import build                                       # noqa: E402
T = float(sys.argv[1]) if len(sys.argv) > 1 else 13.0
P, _ = load_params(); E18 = extraction(P)
P["extraction"]["ve_energy_MeV"] = T; E = extraction(P)
parts, _ = build(P, E, "neutronics")
pv.OFF_SCREEN = True
pl = pv.Plotter(window_size=(1600, 1300)); pl.set_background("white")
col = {"steel_1010": "#9aa0a6", "copper_coil_homog": "#b87333", "copper_ofe": "#d08a50", "stainless_304L": "#c8ccd2",
       "niobium": "#7d86a8", "aluminum_6061": "#e2e4e8"}
tmp = pathlib.Path("/tmp/kiube_ve"); tmp.mkdir(exist_ok=True)
i = E["ve"]["ports"][0]["port"]
tx, ty = E["targets"][i]["x"], E["targets"][i]["y"]
for p in parts:
    bb = p["shape"].BoundingBox()
    if bb.zmin > 1 or bb.zmax < -1 or p["mat"] in ("vacuum",):
        continue
    sec = p["shape"].intersect(cq.Solid.makeBox(4000, 4000, 2, pnt=cq.Vector(-2000, -2000, -1)))
    if sec.Volume() < 1e-3:
        continue
    f = tmp / f"{p['name']}.stl"; cq.exporters.export(sec, str(f), tolerance=0.2)
    pl.add_mesh(pv.read(f), color=("#3b6fd8" if "wedge" in p["name"] else "#e05a2b" if "stripper" in p["name"]
                                   else col.get(p["mat"], "#bbbbbb")), show_edges=False)
for path, c, lab in ((E18["paths"][i], "#1a9a3a", "18 MeV"), (E["paths"][i], "#d0182a", f"{T:g} MeV")):
    pts = np.c_[path, np.full(len(path), 2.0)]
    pl.add_mesh(pv.lines_from_points(pts), color=c, line_width=5, label=lab)
pl.add_legend(bcolor="white", size=(0.16, 0.1), loc="lower left")
pl.add_text(f"Energía variable {T:g} MeV — puerto #8: horquilla larga (naranja) y brida en cuña (azul) {E['ve']['ports'][0]['wedge_deg']:.1f}°",
            position="upper_left", font_size=12, color="black")
pl.view_xy(); cxf, cyf = 0.62 * tx, 0.75 * ty; pl.camera.focal_point = (cxf, cyf, 0); pl.camera.position = (cxf, cyf, 3000)
pl.camera.parallel_projection = True; pl.camera.parallel_scale = 330
out = ROOT / "validation" / f"energia_variable_{T:g}MeV.png"
pl.screenshot(str(out)); print(out, E["ve"]["ports"][0])
