"""Benchmark de blindaje frente al SRG de IBA (MID 55672 Rev. G, Tablas 8 y 9).
Genera: kiube_srg.gdml (ciclotrón, LOD neutrónico, sin búnker; mundo amplio) y srg_config.txt para srg_bench (Geant4).
El búnker estándar del SRG (4 × 4 × 2.5 m interior, hormigón 2.35 g/cm³) se construye en C++ con espesor variable."""
import json, math, pathlib, sys, copy
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import yaml                                                                        # noqa: E402
from kiube import load_params                                                      # noqa: E402
from kiube.beam import extraction                                                  # noqa: E402
from kiube.mc import mc_volumes, conformal_surface_mesh, write_gdml, target_sources  # noqa: E402

PUERTAS = "--puertas" in sys.argv                          # v0.27: el SRG (Tablas 8-9) es sin puertas; con --puertas se añaden las de HDPE
OUT = ROOT / ("out/mc_srg_puertas" if PUERTAS else "out/mc_srg"); OUT.mkdir(parents=True, exist_ok=True)
P, raw = load_params(); E = extraction(P)
P = copy.deepcopy(P); P["vault"]["enabled"] = False
if not PUERTAS:
    P["target"]["puertas_blindadas"] = []
mats = yaml.safe_load(open(ROOT / "materials.yaml"))
vols, info = mc_volumes(P, E, "openmc")
coords, meshes = conformal_surface_mesh(vols)
g = write_gdml(vols, coords, meshes, mats, OUT / "kiube_srg.gdml", world_margin=3700.0)
S = target_sources(P, E)                                    # 180 µA, doble haz (PROP); posiciones con la cavidad real
spec = json.load(open(ROOT / "data/h2o18_spectrum_bon21.json"))
H = P["yoke"]["height"] / 2
floor_z = H - P["yoke"]["floor_to_top"]
c = mats["concrete_portland"]
with open(OUT / "srg_config.txt", "w") as f:
    f.write(f"floor_z {floor_z}\n")
    f.write(f"vault_inner {P['sitio']['vault_min'][0] + 200} {P['sitio']['vault_min'][1] + 200} {P['sitio']['vault_min'][2]}\n")  # 4000×4000×2500 (recomendado)
    f.write(f"concrete_density {c['density']}\n")
    f.write("concrete_wt " + " ".join(f"{k} {v}" for k, v in c["wt"].items()) + "\n")
    for s in S["sources"]:
        if s["type"] == "F18":
            a = math.radians(s["beam_dir_deg"])
            f.write(f"source {s['name']} {s['pos_mm'][0]:.3f} {s['pos_mm'][1]:.3f} {s['pos_mm'][2]:.3f} {math.cos(a):.6f} {math.sin(a):.6f}\n")
    Es, w = spec["E_MeV"], spec["EdNdE_per_proton"]
    o = sorted(range(len(Es)), key=lambda i: Es[i])
    f.write(f"spectrum {len(o)}\n")
    for i in o:
        f.write(f"{Es[i]:.6e} {w[i]:.6e}\n")
json.dump(dict(gdml=g, sources=[s for s in S["sources"] if s["type"] == "F18"], floor_z=floor_z,
               yield_BON21=spec["yield_n_per_proton"], yield_SRG=P["fuente_rad"]["n_por_p_18O"]),
          open(OUT / "geometria.json", "w"), indent=1, default=float)
print(g, floor_z, [(s["name"], s["pos_mm"], s["beam_dir_deg"]) for s in S["sources"] if s["type"] == "F18"])
