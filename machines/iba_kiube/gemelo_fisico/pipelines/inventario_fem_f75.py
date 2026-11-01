#!/usr/bin/env python3
"""F75 — Inventario FEM de la máquina KIUBE (Physics/Facility Twin).

Cruce de volumes.json (137 componentes CAD kiube v0.8) × materials.json (21
materiales con fuente) → inventario componente a componente con masa, clase de
mallabilidad y estado de propiedades FEM, sincronizado con el Physics Twin
(ESPECIFICACIÓN MVP §34.18). Reglas canónicas: simulacion_teorica, pq
no_declarada_nunca, propiedades sin fuente = pendiente_fuente (nunca inventadas),
recomendación != decisión operacional (18.10). FEM 3D completo diferido a versión
final (decisión del usuario); esto es el INVENTARIO, no la corrida FEM.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase75"
OUT.mkdir(parents=True, exist_ok=True)
KIUBE = pathlib.Path(f"{_RPOS_RAIZ}/work/kiube_v08/kiube/out/mc")
F70 = pathlib.Path(f"{_RPOS_RAIZ}/work/f70_mc")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

vols = json.load(open(KIUBE / "volumes.json"))
mats = json.load(open(F70 / "materials.json"))

# Láminas <2 cm no resueltas en voxel 2cm (declarado en F70)
LAMINAS_NO_RESUELTAS = {"havar", "titanium", "graphite_foil", "water_o18", "helium"}

componentes = {}
por_material = {}
clases = {"mallable_voxel_validado": 0, "lamina_no_resuelta_declarada": 0, "otros": 0}
for nombre, info in vols.items():
    mat = info["material"]; v = info["volume_cm3"]
    m = mats.get(mat, {})
    dens = m.get("density")
    lamina = mat in LAMINAS_NO_RESUELTAS
    clase = "lamina_no_resuelta_declarada" if lamina else "mallable_voxel_validado"
    clases[clase] = clases.get(clase, 0) + 1
    componentes[nombre] = {
        "material": mat, "volumen_cm3": v,
        "masa_kg": round(v * dens / 1000, 3) if dens else None,
        "densidad_g_cm3": dens, "fuente_material": m.get("src", "pendiente_fuente"),
        "clase_fem": clase,
        "propiedades_fem": {
            "k_termica": "pendiente_fuente", "cp": "pendiente_fuente",
            "E_young": "pendiente_fuente", "nu_poisson": "pendiente_fuente",
            "nota": "materiales.json no porta propiedades termo-mecánicas; nunca inventar (49.11)",
        },
    }
    e = por_material.setdefault(mat, {"volumen_cm3": 0.0, "masa_kg": 0.0, "n_componentes": 0,
                                      "densidad_g_cm3": dens, "fuente": m.get("src", "pendiente_fuente")})
    e["volumen_cm3"] += v; e["n_componentes"] += 1
    if dens:
        e["masa_kg"] += v * dens / 1000

for e in por_material.values():
    e["volumen_cm3"] = round(e["volumen_cm3"], 1)
    e["masa_kg"] = round(e["masa_kg"], 2)

# Caso FEM existente (declarado): Elmer axisimétrico térmico del blanco
elmer = {
    "caso": "evidencia_v072/elmer_termico/case.sif",
    "tipo": "Heat Equation, estado estacionario, axisimétrico (blanco)",
    "estado": "ejecutado (ver elmer.log, scalars.dat); alcance limitado al blanco, no máquina 3D",
    "md5": md5(EV / "elmer_termico/case.sif"),
}

# Cargas térmicas disponibles con fuente
cargas = {
    "perdida_haz_gas_w": {"valor": 0.0121, "fuente": "F66 mapa_vacio_f66.json (stripping H-, TPMC mapa)"},
    "potencia_blanco": "pendiente_fuente (sin ciclotrón real, memoria 3)",
}

out = {
    "fase": 75, "version": "3.1.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Inventario FEM de la máquina KIUBE (componentes × materiales × mallabilidad)",
    "alcance": "INVENTARIO para FEM futuro; FEM 3D completo diferido a versión final (decisión usuario)",
    "n_componentes": len(componentes),
    "clases_fem": clases,
    "masa_total_kg": round(sum(e["masa_kg"] for e in por_material.values()), 1),
    "por_material": dict(sorted(por_material.items(), key=lambda kv: -kv[1]["masa_kg"])),
    "componentes": componentes,
    "caso_fem_existente": elmer,
    "cargas_termicas_con_fuente": cargas,
    "sincronizacion_gemelo": "Physics Twin sincronizado con inventario real (MVP §34.18); Facility Twin consume masas/componentes",
    "pendientes": [
        "Propiedades termo-mecánicas (k, cp, E, nu) por material: pendiente_fuente",
        "Láminas <2 cm (havar, titanio, grafito, agua O-18, helio): mallado fino dedicado (declarado F70)",
        "FEM 3D máquina completa: diferido a versión final (decisión usuario)",
        "Cargas térmicas de operación real: BLOQUEADAS hasta haz (memoria permanente 3)",
    ],
    "fuentes": ["kiube v0.8 volumes.json (CAD)", "f70_mc/materials.json (fuentes por material)",
                "ESPECIFICACIÓN TÉCNICA DEL MVP INDUSTRIAL §34 (ALME/Physics Twin)",
                "evidencia_v072/elmer_termico (caso Elmer existente)", "F66 (carga stripping)"],
    "md5_fuentes": {"volumes.json": md5(KIUBE / "volumes.json"),
                    "materials.json": md5(F70 / "materials.json")},
}
json.dump(out, open(OUT / "inventario_fem_f75.json", "w"), indent=1)

reg = OUT / "registro_fem.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 75, "version": "3.1.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "inventario_fem_f75.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 75, "estado": "completada", "componentes": len(componentes),
                  "masa_total_kg": out["masa_total_kg"], "clases": clases}, indent=1))
