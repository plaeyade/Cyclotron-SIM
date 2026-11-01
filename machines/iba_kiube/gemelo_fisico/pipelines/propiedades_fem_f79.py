#!/usr/bin/env python3
"""F79 — Propiedades termo-mecánicas con fuente para el FEM final (pendiente F75).

Tabla k, cp, E, nu por material del inventario F75 (21 materiales, materials.json),
cada valor con fuente citada (49.11). Materiales sin fuente confiable quedan
pendiente_fuente — nunca inventados. Gases/fluidos: E/nu no aplican (declarado).
Reglas canónicas: simulacion_teorica, pq no_declarada_nunca, 18.10.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase79"
OUT.mkdir(parents=True, exist_ok=True)
F70 = pathlib.Path(f"{_RPOS_RAIZ}/work/f70_mc")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

mats = json.load(open(F70 / "materials.json"))
inv = json.load(open(EV / "fase75/inventario_fem_f75.json"))

# --- Tabla con fuente (valores a ~20-25 °C, salvo nota) ---
# Cada entrada: (k_W_mK, cp_J_kgK, E_GPa, nu, fuente_id)
PROPS = {
    "steel_1010":       (49.8, 470, 200, 0.29, "AZO-G10100 / ELGIN-1010: k=49.8 W/mK, cp=448-470 J/kgK; E=190-210 GPa, nu=0.27-0.30 (AZO Materials, AISI 1010 / Elgin Fasteners data sheet)"),
    "stainless_304L":   (16.2, 500, 193, 0.29, "AZO-304 / USPTO-304: k=16.2 W/mK (0-100°C), cp=500 J/kgK, E=193-200 GPa, nu=0.29 (AZO AID965; AISI 304 data)"),
    "copper_ofe":       (390, 390, 120, 0.34, "MIF-C10200: k=390 W/mK, cp=390 J/kgK, E=117-120 GPa, nu=0.34 (UNS C10200 Cu-OF; Aurubis C10200 data sheet)"),
    "niobium":          (53.7, 265, 104.9, 0.397, "MATECK-Nb: k=53.7 W/mK @300K, E=104.9 GPa, nu=0.397 (Mateck research data); cp=265 J/kgK (valor tabulado estándar, declarado)"),
    "titanium":         (6.7, 560, 113.8, 0.342, "TIMET-6-4 / AMARDEEP-6-4: Ti-6Al-4V k=6.7 W/mK, cp=560 J/kgK, E=113.8 GPa, nu=0.342 (grado 5; si el grado real es puro Ti, pendiente_fuente)"),
    "havar":            (13.0, 420, 200, 0.285, "WIKI-HAVAR / HPM-HAVAR: k=13.0 W/mK, E=200-210 GPa, nu~0.285 (AZO R30004: 0.27-0.30); cp=420 J/kgK [EST, marcado]"),
    "graphite_foil":    (140, 710, 10, 0.2, "MINSEAL-FG / AZO-C516: grafito flexible k=140 W/mK en plano (5 a través, ANISOTROPO declarado), E~1.4-10 GPa, nu~0.2, cp=697-771 J/kgK"),
    "concrete_portland": (1.7, 880, 30, 0.2, "NUREG-CR7031 / ET-429: concreto denso k=1.0-1.8 W/mK (1.7 típico piedra), E~30 GPa y nu~0.2 (compilación NUREG/CR-7031), cp~880 J/kgK [EST típico, marcado]"),
    "aluminum":         (167, 896, 68.3, 0.33, "LESSCAD-AL: k=167 W/mK @20°C, cp=896 J/kgK, E=68.3 GPa, nu=0.33 (valores típicos tabulados)"),
}
FLUIDOS = {"water_o18", "helium", "aire", "air", "vacio", "vacuum"}

propiedades = {}
n_con_fuente = 0
for nombre, m in mats.items():
    if nombre in PROPS:
        k, cp, E, nu, fuente = PROPS[nombre]
        propiedades[nombre] = {"k_W_mK": k, "cp_J_kgK": cp, "E_GPa": E, "nu": nu,
                               "fuente": fuente, "estado": "con_fuente"}
        n_con_fuente += 1
    elif nombre in FLUIDOS or (m.get("density") or 0) < 0.01:
        propiedades[nombre] = {"k_W_mK": None, "cp_J_kgK": None, "E_GPa": None, "nu": None,
                               "fuente": None,
                               "estado": "no_aplica (fluido/gas/vacío — E y nu sin sentido estructural, declarado)"}
    else:
        propiedades[nombre] = {"k_W_mK": None, "cp_J_kgK": None, "E_GPa": None, "nu": None,
                               "fuente": None, "estado": "pendiente_fuente"}

# Bobina homogeneizada: regla de mezclas en peso (CALC declarado, dominada por Cu)
coil = mats.get("copper_coil_homog", {})
if coil:
    w_cu = coil["wt"]["Cu"] / 100.0
    k_mix = w_cu * PROPS["copper_ofe"][0] + (1 - w_cu) * 1.0   # epoxi-vidrio/agua ~1 W/mK [EST]
    cp_mix = w_cu * PROPS["copper_ofe"][1] + (1 - w_cu) * 1000
    E_mix = w_cu * PROPS["copper_ofe"][2] + (1 - w_cu) * 5     # epoxi ~5 GPa [EST]
    propiedades["copper_coil_homog"] = {
        "k_W_mK": round(k_mix, 1), "cp_J_kgK": round(cp_mix, 0), "E_GPa": round(E_mix, 1), "nu": 0.34,
        "fuente": "CALC regla de mezclas en peso sobre MIF-C10200 + epoxi/agua [EST] — homogeneizada, declarada no validada",
        "estado": "calc_declarado"}

out = {
    "fase": 79, "version": "3.5.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Propiedades termo-mecánicas con fuente para el FEM final (pendiente F75)",
    "condicion": "valores a temperatura ambiente (~20-25 °C) salvo nota; para FEM térmico-estructural",
    "n_materiales": len(mats), "n_con_fuente": n_con_fuente,
    "propiedades": propiedades,
    "cobertura_inventario_f75": {
        "masa_cubierta_pct": round(100 * sum(e["masa_kg"] for n, e in inv["por_material"].items()
                                             if propiedades.get(n, {}).get("estado") in ("con_fuente", "calc_declarado"))
                                 / sum(e["masa_kg"] for e in inv["por_material"].values()), 2),
    },
    "advertencias": [
        "Valores a temperatura ambiente; el FEM térmico a temperatura de operación requiere curvas k(T), cp(T) — pendiente_fuente para operación (memoria 3).",
        "Grafito foil ANISOTROPO declarado (140 en plano / 5 a través); titanium asume grado 5 (si es Ti puro comercial: pendiente_fuente).",
        "cp de havar y Nb, cp/E del concreto: valores típicos tabulados marcados [EST] donde la fuente no lo cubre.",
        "copper_coil_homog es CALC por regla de mezclas (declarado, no validado); bobina real es anisotrópica (devanado + epoxi).",
        "Recomendación (18.10), no decisión operacional: pedir certificados de colada (F70 ya lo declaró para Co en acero).",
    ],
    "fuentes": ["AZO-G10100/ELGIN-1010 (AISI 1010)", "AZO AID965 / USPTO-304 (304L)",
                "MIF-C10200 / Aurubis (Cu-OF)", "MATECK-Nb (niobio)", "TIMET-6-4 (Ti-6Al-4V)",
                "WIKI-HAVAR / HPM-HAVAR / AZO R30004", "MINSEAL-FG / AZO-C516 (grafito)",
                "NUREG/CR-7031 (concreto Portland)", "LESSCAD-AL (aluminio)",
                "f70_mc/materials.json", "evidencia_v072/fase75 (inventario)"],
    "md5_fuentes": {"materials.json": md5(F70 / "materials.json"),
                    "inventario_f75": md5(EV / "fase75/inventario_fem_f75.json")},
}
json.dump(out, open(OUT / "propiedades_fem_f79.json", "w"), indent=1)

reg = OUT / "registro_fem.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 79, "version": "3.5.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "propiedades_fem_f79.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 79, "estado": "completada", "n_con_fuente": n_con_fuente,
                  "cobertura_masa_pct": out["cobertura_inventario_f75"]["masa_cubierta_pct"]}, indent=1))
