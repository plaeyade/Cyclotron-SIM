"""Exporta los dominios del sector de 90° (fem_mag) con parámetros modificados, para estudios de
sensibilidad con fem/magnet_ngsolve.py --dom. Uso:
  python fem/fem_mag_variante.py --tag valle100 --set pole.valley_floor_z=100"""
import argparse, json, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import cadquery as cq
from kiube import load_params
from kiube.beam import extraction
from kiube.geometry import build_fem_mag
ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--set", nargs="*", default=[])
a = ap.parse_args()
P, _ = load_params()
for kv in a.set:
    k, v = kv.split("="); g, kk = k.split(".")
    if v.startswith("["):
        P[g][kk] = json.loads(v)
        continue
    try:
        P[g][kk] = float(v)
    except ValueError:
        P[g][kk] = v
E = extraction(P)
d = ROOT / "out/fem_mag_variantes" / a.tag / "sector90"
(d / "domains").mkdir(parents=True, exist_ok=True)
doms = build_fem_mag(P, E, sector=True)
for p in doms:
    cq.exporters.export(p["shape"], str(d / "domains" / f"{p['name']}.step"))
json.dump([{"name": p["name"], "material": p["mat"], "role": p["role"], "volume_mm3": round(p["shape"].Volume(), 1)}
           for p in doms], open(d / "domains.json", "w"), indent=1)
json.dump({"set": a.set}, open(d.parent / "variante.json", "w"))
print("ok", d)
