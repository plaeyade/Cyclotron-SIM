"""Configuración del equivalente v0.35 de RPOS F70 (mapa de dosis y activación de la máquina en el búnker del proyecto)
para mc_benchmark/bunker_mapa. Bloque de hormigón igual al de F70 para poder comparar: muros 2.11 m, losa de piso 0.55 m,
techo 1.45 m (F70: ±506 × ±446 cm, z −167…283 cm); malla H*(10) de 10 cm que cubre el bloque + 1 m.
(Cabecera original de generar_techo.py:) Configuración del cálculo del TECHO del búnker del proyecto (v0.28).

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
ap.add_argument("--losa", type=float, default=550.0, help="losa de piso [mm] (F70: 550)")
ap.add_argument("--techo", type=float, default=1450.0, help="techo [mm] (F70: 1450; SRG mín. 1700)")
a = ap.parse_args()
P, raw = load_params(); E = extraction(P)
P = copy.deepcopy(P); P["vault"]["enabled"] = False
if a.sin_puertas:
    P["target"]["puertas_blindadas"] = []
V = {k: (v["v"] if isinstance(v, dict) and "v" in v else v) for k, v in P["vault"].items()}
ih = a.altura or V["inner_height"]
OUT = ROOT / ("out/mc_f70v035_sinpuertas" if a.sin_puertas else "out/mc_f70v035") / f"h{int(ih)}_t{int(a.techo)}"
OUT.mkdir(parents=True, exist_ok=True)
mats = yaml.safe_load(open(ROOT / "materials.yaml"))
vols, info = mc_volumes(P, E, "openmc")
coords, meshes = conformal_surface_mesh(vols)
g = write_gdml(vols, coords, meshes, mats, OUT / "kiube.gdml", world_margin=6500.0)   # cabe el bloque + 1 m de malla
S = target_sources(P, E)
spec = json.load(open(ROOT / "data/h2o18_spectrum_bon21.json"))
H = P["yoke"]["height"] / 2
floor_z = H - P["yoke"]["floor_to_top"]
c = mats["concrete_portland"]
with open(OUT / "mapa_config.txt", "w") as f:
    f.write(f"floor_z {floor_z}\n")
    f.write(f"vault_inner {V['inner_x']} {V['inner_y']} {ih}\n")
    f.write(f"wall {V['wall']}\n")
    f.write(f"slab {a.losa}\nroof {a.techo}\n")
    r100 = lambda v: math.ceil(v / 100.0) * 100.0                       # malla de 10 cm alineada
    hx, hy = r100(V["inner_x"] / 2 + V["wall"] + 1000), r100(V["inner_y"] / 2 + V["wall"] + 1000)
    z0, z1 = -r100(-(floor_z - a.losa - 1000)), r100(floor_z + ih + a.techo + 1000)
    f.write(f"mesh {-hx} {hx} {-hy} {hy} {z0} {z1} 100\n")
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
               wall=V["wall"], losa=a.losa, techo=a.techo, puertas=not a.sin_puertas, version_modelo="v0.35",
               yield_BON21=spec["yield_n_per_proton"], yield_SRG=P["fuente_rad"]["n_por_p_18O"]),
          open(OUT / "geometria.json", "w"), indent=1, default=float)
print(OUT, floor_z, [(s["name"], [round(x) for x in s["pos_mm"]]) for s in S["sources"] if s["type"] == "F18"])
