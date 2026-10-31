"""Prepara la escena 3D de la región central para Blender: electrodos (STEP → STL) y órbitas calculadas.
Uso: python fem/region_central_escena.py out/central_region/c18_9
Escribe <dir>/escena/{*.stl, escena.json}; luego:
  blender -b -P blender/render_region_central.py -- --escena <dir>/escena/escena.json --salida rc.png [--gpu auto]"""
import json, pathlib, sys
import numpy as np
import cadquery as cq

d = pathlib.Path(sys.argv[1]); R = json.load(open(d / "resultado.json"))
out = d / "escena"; out.mkdir(exist_ok=True)
MAT = {"colinas_tapon": "acero", "dummy_dees": "acero_oscuro", "chimeneas": "cuw", "dees": "cobre", "pullers": "cobre"}
piezas = []
for grp, files in R.get("electrodos_step", {}).items():
    for f in files:
        try:
            sh = cq.importers.importStep(str(d / "electrodos" / f))
        except Exception as e:                       # noqa: BLE001
            print("[aviso]", f, e); continue
        stl = out / (pathlib.Path(f).stem + ".stl")
        cq.exporters.export(sh, str(stl), tolerance=0.05, angularTolerance=0.1)
        piezas.append(dict(stl=stl.name, material=MAT.get(grp, "acero"), grupo=grp))
Z = np.load(d / "trayectorias.npz")
orb = []
COLOR = {"H-": [0.25, 1.0, 0.15], "D-": [0.10, 0.22, 1.0]}
for e in R["especies"]:
    key = e["especie"].replace("-", "m")
    dz = e["fase_diseno"]
    if dz:
        t = Z[f"{key}_fase_{dz['fase']:+.1f}"]
        orb.append(dict(especie=e["especie"], tipo="diseno", color=COLOR[e["especie"]], xyz=t[::2, 1:4].round(3).tolist()))
    for k in range(40):
        n = f"{key}_haz_{k:03d}"
        if n in Z.files and len(Z[n]) > 50:
            orb.append(dict(especie=e["especie"], tipo="haz", color=COLOR[e["especie"]], xyz=Z[n][::4, 1:4].round(3).tolist()))
    mk = e.get("marcadores_RF") or {}
    for tipo, col in (("V0", [0.9, 0.08, 0.08]), ("Vmax", [0.1, 0.8, 0.2])):
        for p in mk.get(tipo, []):
            orb.append(dict(especie=e["especie"], tipo="marcador", color=col, xyz=[[p[0], p[1], 0.0]]))
json.dump(dict(piezas=piezas, orbitas=orb, unidades="mm"), open(out / "escena.json", "w"))
print(len(piezas), "electrodos,", len(orb), "curvas/marcadores →", out / "escena.json")
