"""Configuración del cálculo del TECHO del búnker del proyecto (v0.28) para mc_benchmark/bunker_techo.

Búnker del proyecto (RPOS, plano BIND 02-10-26): interior 5.9 × 4.7 m (estimado sobre planta), muros 2.11 m;
altura interior y espesor del techo pendientes de los cortes → se barren. Encima: cuarto de aire acondicionado.
Fuente: blancos F18 T1 y T5 en la cara del yugo, CON las puertas de HDPE (configuración del proyecto); --sin_puertas
las quita. El ciclotrón se centra en la sala (posición real pendiente de los cortes).
Uso: python mc_benchmark/generar_techo.py [--sin_puertas] [--altura 2500]  → out/mc_techo[_sinpuertas]/
"""
import argparse, copy, json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import yaml                                                                        # noqa: E402
from kiube import load_params                                                      # noqa: E402
from kiube.beam import extraction                                                  # noqa: E402
from kiube.mc import mc_volumes, conformal_surface_mesh, write_gdml, target_sources  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--sin_puertas", action="store_true")
ap.add_argument("--altura", type=float, default=None, help="altura interior del búnker [mm] (defecto: vault.inner_height)")
a = ap.parse_args()
P, raw = load_params(); E = extraction(P)
P = copy.deepcopy(P); P["vault"]["enabled"] = False
if a.sin_puertas:
    P["target"]["puertas_blindadas"] = []
V = {k: (v["v"] if isinstance(v, dict) and "v" in v else v) for k, v in P["vault"].items()}
ih = a.altura or V["inner_height"]
OUT = ROOT / ("out/mc_techo_sinpuertas" if a.sin_puertas else "out/mc_techo") / f"h{int(ih)}"
OUT.mkdir(parents=True, exist_ok=True)
mats = yaml.safe_load(open(ROOT / "materials.yaml"))
vols, info = mc_volumes(P, E, "openmc")
coords, meshes = conformal_surface_mesh(vols)
g = write_gdml(vols, coords, meshes, mats, OUT / "kiube.gdml", world_margin=6500.0)   # cabe muro 2.11 m + techo + detector
S = target_sources(P, E)
spec = json.load(open(ROOT / "data/h2o18_spectrum_bon21.json"))
H = P["yoke"]["height"] / 2
floor_z = H - P["yoke"]["floor_to_top"]
c = mats["concrete_portland"]
with open(OUT / "techo_config.txt", "w") as f:
    f.write(f"floor_z {floor_z}\n")
    f.write(f"vault_inner {V['inner_x']} {V['inner_y']} {ih}\n")
    f.write(f"wall {V['wall']}\n")
    f.write(f"concrete_density {c['density']}\n")
    f.write("concrete_wt " + " ".join(f"{k} {v}" for k, v in c["wt"].items()) + "\n")
    for s in S["sources"]:
        if s["type"] == "F18":
            b = math.radians(s["beam_dir_deg"])
            f.write(f"source {s['name']} {s['pos_mm'][0]:.3f} {s['pos_mm'][1]:.3f} {s['pos_mm'][2]:.3f} {math.cos(b):.6f} {math.sin(b):.6f}\n")
    Es, w = spec["E_MeV"], spec["EdNdE_per_proton"]
    o = sorted(range(len(Es)), key=lambda i: Es[i])
    f.write(f"spectrum {len(o)}\n")
    for i in o:
        f.write(f"{Es[i]:.6e} {w[i]:.6e}\n")
json.dump(dict(gdml=g, sources=[s for s in S["sources"] if s["type"] == "F18"], floor_z=floor_z, inner=[V["inner_x"], V["inner_y"], ih],
               wall=V["wall"], puertas=not a.sin_puertas,
               yield_BON21=spec["yield_n_per_proton"], yield_SRG=P["fuente_rad"]["n_por_p_18O"]),
          open(OUT / "geometria.json", "w"), indent=1, default=float)
print(OUT, floor_z, [(s["name"], [round(x) for x in s["pos_mm"]]) for s in S["sources"] if s["type"] == "F18"])
