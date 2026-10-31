"""Barrido de RF v0.19: efecto del puente central y del sintonizador (casos en out/fem_rf/resultado_<tag>.json).
Uso: python fem/rf_barrido_v19.py tag clave=valor [clave=valor ...]"""
import copy, json, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
from rf_kiube import resolver
from kiube import load_params
from kiube.beam import extraction
tag = sys.argv[1]
P, _ = load_params(); P = copy.deepcopy(P)
for kv in sys.argv[2:]:
    k, v = kv.split("="); g, n = k.split(".")
    P[g][n] = float(v)
E = extraction(P)
import os
res, _ = resolver(P, E, maxh=float(os.environ.get("MAXH", 0.045)), grading=float(os.environ.get("GRAD", 0.55)), orden=2, tag=tag, n_modos=4, shift_frac=0.5)
res.pop("campo_en_fase", None); res.pop("escala_40kV", None)
res["param"] = sys.argv[2:]
json.dump(res, open(ROOT / f"out/fem_rf/resultado_{tag}.json", "w"), indent=1, default=str)
for m in res["modos"][:2]:
    print(f"RES {tag} {m['fase']} {m['f_MHz']:.3f} MHz V1/V0 {m['V1_sobre_V0']:.3f} Q {m['Q_cobre']:.0f} P {m['P_kW_40kV']:.2f} kW")
