"""Modelo frente al Site Requirement Guide de IBA (MID 55672 Rev. G, 2023) y la cotización del proyecto (PROP).
Comprueba cotas de instalación (altura total, plano medio, pies, envolvente) y la masa total.
Uso: python validation/validar_srg.py   → validation/validacion_srg.json
NOTA: el SRG es un documento propietario de IBA; aquí solo se usan cifras de instalación para verificar el modelo."""
import json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import yaml                                                                     # noqa: E402
from kiube import load_params                                                   # noqa: E402
from kiube.beam import extraction                                               # noqa: E402
from kiube.geometry import build                                                # noqa: E402

P, raw = load_params(); E = extraction(P)
mats = yaml.safe_load(open(ROOT / "materials.yaml"))
H = P["yoke"]["height"] / 2
floor = H - P["yoke"]["floor_to_top"]
res = {"documentos": {"SRG": "MID 55672 Rev. G (2023), propietario IBA",
                      "PROP": "cotización IBA 2025/8499 v02.01 (solo configuración técnica)"}, "checks": []}


def chk(nombre, modelo, ref, tol, unidad, fuente, nota=""):
    ok = abs(modelo - ref) <= tol
    res["checks"].append(dict(nombre=nombre, modelo=round(modelo, 1), referencia=ref, tolerancia=tol, unidad=unidad,
                              ok=bool(ok), fuente=fuente, nota=nota))
    print(f"[{'OK' if ok else '!!'}] {nombre}: modelo {modelo:.1f} vs {ref} ±{tol} {unidad}  ({fuente})")


parts, _ = build(P, E, "render")
bb = None
for p in parts:
    if p["name"].startswith("self_shield"):
        continue
    b = p["shape"].BoundingBox()
    bb = [min(bb[0], b.xmin), min(bb[1], b.ymin), min(bb[2], b.zmin), max(bb[3], b.xmax), max(bb[4], b.ymax),
          max(bb[5], b.zmax)] if bb else [b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax]
chk("plano medio sobre el piso", -floor, 1120.0, 5.0, "mm", "SRG §4.3.5 / SRG-0420")
chk("altura total (cima del cáncamo)", bb[5] - floor, 1850.0, 5.0, "mm", "SRG Tabla 1 / SRG-0080")
chk("piso del modelo (pie más bajo)", bb[2] - floor, 0.0, 1.0, "mm", "SRG")
import re as _re                                              # v0.27: los blancos (T1…T8) van en la cara del yugo y se montan
yoke_parts = [p for p in parts if not p["name"].startswith(("self_shield", "elec_", "hose_", "pneumatic", "water_h",   # en sitio:
                                                            "water_b", "water_v", "flow_tx"))      # no cuentan para la
              and not _re.match(r"^T\d_", p["name"])]                                                 # envolvente de maniobra
bx = None
for p in yoke_parts:
    b = p["shape"].BoundingBox()
    bx = [min(bx[0], b.xmin), min(bx[1], b.ymin), max(bx[2], b.xmax), max(bx[3], b.ymax)] if bx else \
        [b.xmin, b.ymin, b.xmax, b.ymax]
for eje, L in (("X", bx[2] - bx[0]), ("Y", bx[3] - bx[1])):    # entre la ficha (1900) y el máximo de maniobra (2050)
    chk(f"envolvente {eje} (sin servicios ni blancos)", L, 1975.0, 75.0, "mm",
        "SRG Tabla 1 (2000) y SRG-0080 (máx. 2050); ficha DS 1900")
legs = [p for p in parts if p["name"].startswith("leg_")]
c = legs[0]["shape"].Center()
chk("posición de los pies (|x|, |y|)", (abs(c.x) + abs(c.y)) / 2, 525.0, 1.0, "mm", "SRG-0315")
pb = legs[0]["shape"].BoundingBox()
chk("pie en el piso (lado)", pb.xmax - pb.xmin, 350.0, 1.0, "mm", "SRG-0315")

# masa: piezas reales (fem_struct = hardware ≥ 0.5 mm) + bobina; sin cosméticos ni fluidos
pm, _ = build(P, E, "fem_struct")
m = sum(p["shape"].Volume() * 1e-6 * mats[p["mat"]]["density"] for p in pm)
chk("masa del ciclotrón", m, 18000.0, 1800.0, "kg", "SRG Tabla 1/4 (18 000 kg)",
    "tolerancia ±10 %: faltan piezas internas (RF, bombas reales, cableado) y hay cotas EST")
res["masa_por_material_kg"] = {}
for p in pm:
    res["masa_por_material_kg"][p["mat"]] = round(res["masa_por_material_kg"].get(p["mat"], 0) +
                                                  p["shape"].Volume() * 1e-6 * mats[p["mat"]]["density"], 1)
# configuración del proyecto
TG = P["target"]
res["proyecto"] = dict(modelo=P["proyecto"]["modelo"], blancos=TG["types"],
                       strippers=P["extraction"]["strippers_instalados"],
                       strippers_modelo=sorted({int(p["name"].split("_")[-1]) for p in parts
                                                if p["name"].startswith("stripper_arm_")}))
ok = res["proyecto"]["strippers"] == res["proyecto"]["strippers_modelo"]
res["checks"].append(dict(nombre="strippers instalados = PROP (4 de serie)", ok=ok, fuente="PROP 201.052"))
print(f"[{'OK' if ok else '!!'}] strippers instalados {res['proyecto']['strippers_modelo']} (PROP: 4 de serie)")
from kiube.geometry import target_cavity                                         # noqa: E402
res["blancos_volumen_ml"] = {t: round(math.pi * (target_cavity(TG, t)[0] / 2) ** 2 * target_cavity(TG, t)[1] / 1000, 2)
                             for t in sorted(set(TG["types"]) - {"none"})}
print("volúmenes de cavidad (ml):", res["blancos_volumen_ml"])
res["resumen"] = dict(ok=sum(1 for c in res["checks"] if c["ok"]), total=len(res["checks"]))
json.dump(res, open(ROOT / "validation/validacion_srg.json", "w"), indent=1, ensure_ascii=False)
print(res["resumen"])
