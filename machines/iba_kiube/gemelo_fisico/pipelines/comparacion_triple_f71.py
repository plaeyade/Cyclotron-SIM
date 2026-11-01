#!/usr/bin/env python3
"""F71 — Comparación triple de blindaje: kiube-voxel (F70) vs casco DAGMC (F62) vs CSG BIND (F52).

Reglas canónicas: simulacion_teorica, pq no_declarada_nunca, métricas NO homogéneas
declaradas explícitamente (cada fase midió cantidades distintas); figuras/cifras solo
con fuente citada (49.11); recomendaciones no son decisiones operativas (18.10).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase71"
OUT.mkdir(parents=True, exist_ok=True)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

# --- Fuentes ---
f52 = json.load(open(EV / "fase52/combinado_sab_f52.json"))           # CSG BIND, transmisión S(a,b)
f62 = json.load(open(EV / "fase62/corrida_mc_f62.json"))              # casco DAGMC
f70e = json.load(open(EV / "fase70/mc_maquina_bunker_f70.json"))      # kiube voxel + búnker BIND

full = f70e["corridas"].get("full_5x2000x2_malla10cm") or {}
verif = f70e["corridas"]["verificacion_1cm_3000x2"]

comparacion = {
    "F52_csg_bind_sab": {
        "geometria": "CSG analítica búnker BIND (RPOS-PLANO), sin máquina",
        "metrica": "transmision a través del muro con S(a,b) (análogo transmisión)",
        "valor": f52["transmision_media"], "sigma": f52["sigma_semillas"],
        "historias": f52["historias_totales"],
        "fuente": "evidencia_v072/fase52/combinado_sab_f52.json",
    },
    "F62_casco_dagmc": {
        "geometria": "casco DAGMC sellado (bunker_dagmc_cm.h5m), sin máquina detallada",
        "metrica": "fraccion_fuga global + flujo en muro",
        "fraccion_fuga": f62["resultados"]["fraccion_fuga"],
        "fraccion_fuga_sigma": f62["resultados"]["fraccion_fuga_sigma"],
        "flujo_muro_cm2_part": f62["resultados"]["flujo_muro_cm_por_particula"],
        "historias": f62["corrida"]["particles"] * f62["corrida"]["batches"],
        "fuente": "evidencia_v072/fase62/corrida_mc_f62.json",
    },
    "F70_kiube_voxel_bunker": {
        "geometria": "máquina KIUBE voxelizada 2 cm + búnker BIND CSG (muro 2.11 m)",
        "metrica": "dosis integral fuera del bloque relativa (verificación) + tallies full",
        "fuera_bloque_rel": verif["dose_integral_fuera_bloque_rel"],
        "dentro_bloque_rel": verif["dose_integral_dentro_bloque_rel"],
        "full_estado": f70e["corridas"]["estado_full"],
        "full_historias": full.get("historias_total"),
        "fuente": "evidencia_v072/fase70/mc_maquina_bunker_f70.json",
    },
}

# --- Alineación de muros (oeste/este) según planta BIND (RPOS-PLANO) ---
alineacion_muros = {
    "fuente_planta": "RPOS-PLANO (BIND), plantas_bind.dxf",
    "muro_principal_espesor_m": 2.11,
    "sala_m": {"x": [-295, 295], "y": [-235, 235], "z": [-125, 125]},
    "bloque_concreto_m": {"x": [-506, 506], "y": [-446, 446], "z": [-167, 283]},
    "criterio": ("Las tres aproximaciones (F52 CSG, F62 DAGMC casco, F70 voxel+CSG) usan la "
                 "misma planta BIND; muro oeste y este definidos por el bloque ±506/±446 cm. "
                 "F70 hereda la bóveda CSG BIND idéntica a F52 (misma definición de muros), "
                 "por lo que oeste/este quedan alineados por construcción; F62 usa casco "
                 "sellado derivado del mismo plano (declarado en fase62)."),
    "estado": "alineados_por_construccion (F52=F70) / derivado_del_mismo_plano (F62)",
}

advertencias = [
    "Métricas NO homogéneas: transmisión (F52), fracción de fuga (F62) y dosis fuera del bloque (F70) "
    "miden cantidades distintas; la comparación es cualitativa, no una verificación cruzada cuantitativa (declarado).",
    "F52 y F62 no incluyen la máquina (autoblindaje); F70 sí (voxel 2 cm validado vs CAD). "
    "Las diferencias de magnitud son esperables y no indican error.",
    "Cifra fuera_bloque de F70 proviene de la corrida de verificación (3000x2, baja estadística, declarado); "
    "la corrida full (5x2000x2) aporta tallies de activación y mallas, no ese cociente.",
    "Recomendación (no decisión operativa, 18.10): si se requiere comparación cuantitativa estricta, "
    "repetir F52/F62 con la métrica de dosis H*(10) idéntica a F70 en fase futura.",
]

out = {
    "fase": 71, "version": "2.7.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca",
    "fecha": datetime.date.today().isoformat(),
    "titulo": "Comparación triple de blindaje: kiube-voxel (F70) vs casco DAGMC (F62) vs CSG BIND (F52)",
    "comparacion": comparacion,
    "alineacion_muros_oeste_este": alineacion_muros,
    "advertencias": advertencias,
    "fuentes": ["RPOS-PLANO (BIND)", "evidencia_v072/fase52", "evidencia_v072/fase62",
                "evidencia_v072/fase70", "ICRP-74/116", "ENDF/B-VIII.0"],
    "md5_fuentes": {k: md5(f"{_RPOS_RAIZ}/work/" + v["fuente"]) for k, v in comparacion.items()},
}
json.dump(out, open(OUT / "comparacion_triple_f71.json", "w"), indent=1)

reg = OUT / "registro_comparacion.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 71, "version": "2.7.0",
                        "fecha": out["fecha"], "md5_evidencia": md5(OUT / "comparacion_triple_f71.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 71, "estado": "completada",
                  "T_F52": comparacion["F52_csg_bind_sab"]["valor"],
                  "fuga_F62": comparacion["F62_casco_dagmc"]["fraccion_fuga"],
                  "fuera_bloque_F70": comparacion["F70_kiube_voxel_bunker"]["fuera_bloque_rel"]}, indent=1))
