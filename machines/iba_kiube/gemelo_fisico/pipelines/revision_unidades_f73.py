#!/usr/bin/env python3
"""F73 — Revisión de unidades de tallies F70 (hallazgo y corrección).

Hallazgo (test A/B, evidencia: diagnostico_ab_strength.json): este build de
OpenMC 0.16.0 normaliza los tallies a la FUERZA DE FUENTE, no por partícula fuente
(razón medida = 2.308e12, exactamente la fuerza de cada fuente T1/T2).
Consecuencia: los valores de f70_mc/full_voxel.json ya son tasas por segundo a
intensidad nominal (4.616e12 n/s, BON21). F72 multiplicó de nuevo por S_n →
inflado por 4.616e12 (declarado allí como verificación pendiente).

Correcciones aplicadas aquí (sin reescribir evidencias históricas, append-only):
  - tasas de activación correctas = mean_sum (ya s^-1)  → A_sat y A(t) EOB
  - verif_voxel.json: reaction_rates_per_s y max_dose_rate estaban ×4.616e12 demás
  - dosis fuera del bloque (cociente) no se ve afectada.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime, math

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase73"
OUT.mkdir(parents=True, exist_ok=True)
F70 = pathlib.Path(f"{_RPOS_RAIZ}/work/f70_mc")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

diag = json.load(open(OUT / "diagnostico_ab_strength.json"))
full = json.load(open(F70 / "full_voxel.json"))
verif = json.load(open(F70 / "verif_voxel.json"))
S_n = 4616232340295.063  # n/s nominal (BON21)

# --- Verificación cuantitativa del mecanismo ---
razon = diag["razon_esperada"]
a = diag["strength_nominal_2.308e12"]["flux_mesh"]["mean_sum"]
b = diag["strength_unitaria"]["flux_mesh"]["mean_sum"]
razon_medida = a / b
desv_rel = abs(razon_medida - razon) / razon

ANO = 365.25 * 86400
ISOTOPOS = {
    "Co60":   {"T12_s": 5.2714 * ANO},
    "Eu152":  {"T12_s": 13.537 * ANO},
    "Ag110m": {"T12_s": 249.79 * 86400},
    "Cu64":   {"T12_s": 12.701 * 3600},
}
CANALES = {
    "steel_1010:Co59(n,gamma)":         "Co60",
    "copper_coil_homog:Cu63(n,a)":      "Co60",
    "concrete_portland:Eu151(n,gamma)": "Eu152",
    "copper_coil_homog:Ag109(n,gamma)": "Ag110m",
    "copper_coil_homog:Cu63(n,gamma)":  "Cu64",
}
ESCENARIOS_S = {"1h": 3600, "8h": 8 * 3600, "30d": 30 * 86400, "1a": ANO}

act = {}
for canal, iso in CANALES.items():
    t = full["tallies"][canal]
    rate = t["mean_sum"]              # CORREGIDO: ya son reacciones/s a fuente nominal
    rate_sd = t["std_sum"]
    lam = math.log(2) / ISOTOPOS[iso]["T12_s"]
    e = act.setdefault(iso, {"A_sat_Bq": 0.0, "A_sat_Bq_var": 0.0, "canales": [],
                             "T12_s": ISOTOPOS[iso]["T12_s"],
                             "fuente_T12": "ENDF/B-VIII.0 decay / NuDat"})
    e["A_sat_Bq"] += rate
    e["A_sat_Bq_var"] += rate_sd ** 2
    e["canales"].append({"canal": canal, "tasa_s^-1": rate, "tasa_sigma": rate_sd})
for iso, e in act.items():
    e["A_sat_Bq_sigma"] = math.sqrt(e.pop("A_sat_Bq_var"))
    lam = math.log(2) / e["T12_s"]
    e["A_EOB_Bq"] = {k: e["A_sat_Bq"] * (1 - math.exp(-lam * ts)) for k, ts in ESCENARIOS_S.items()}

# --- Corrección de cifras históricas de F70 (verif) ---
f72 = json.load(open(EV / "fase72/activacion_gemelo_f72.json"))
correcciones = {
    "verif_reaction_rates_per_s": {
        "estado": "infladas x4.616e12 en evidencia fase70 (multiplicación doble por S_n)",
        "valores_corregidos_s^-1": {
            k: {"mean": v["mean"] / S_n, "std": v["std"] / S_n}
            for k, v in verif["reaction_rates_per_s"].items()
        },
    },
    "verif_max_dose_rate_Sv_h": {
        "original": verif["max_dose_rate_Sv_per_h_in_mesh"],
        "corregido": verif["max_dose_rate_Sv_per_h_in_mesh"] / S_n,
        "nota": "voxel de máxima dosis (zona blanco); cociente fuera/dentro del bloque no afectado",
    },
    "f72_A_sat_Bq": {
        "estado": "infladas x4.616e12 en evidencia fase72 (verificación de unidades declarada pendiente; resuelta aquí)",
        "valores_corregidos": {iso: f72["isotopos"][iso]["A_sat_Bq"] / S_n for iso in f72["isotopos"]},
    },
}

out = {
    "fase": 73, "version": "2.9.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Revisión de unidades de tallies F70: OpenMC 0.16 normaliza a fuerza de fuente",
    "hallazgo": {
        "mecanismo": "tallies normalizados a source strength (no por partícula fuente) en este build de OpenMC 0.16.0",
        "prueba": "test A/B 500x2: misma geometría, strength 2.308e12 vs 1.0",
        "razon_medida_flux_mesh": razon_medida,
        "razon_esperada_strength": razon,
        "desviacion_relativa": desv_rel,
        "evidencia": "evidencia_v072/fase73/diagnostico_ab_strength.json",
    },
    "consecuencia": "full_voxel.json ya contiene tasas por SEGUNDO a fuente nominal 4.616e12 n/s; no multiplicar de nuevo",
    "activacion_corregida": act,
    "dosis_agregada_corregida": {
        "dose_mesh_mean_sum_Sv_cm3_h": full["tallies"]["dose_mesh"]["mean_sum"],
        "nota": "suma de H*(10) sobre malla 10 cm a fuente nominal (ya incluye S_n)",
        "fuera_bloque_rel_verif": verif["dose_integral_fuera_bloque_rel"],
    },
    "correcciones_historicas": correcciones,
    "politica": "evidencias fase70/fase72 NO reescritas (append-only); esta fase suple los valores corregidos",
    "advertencias": [
        "Escenario teórico continuo a fuente nominal; sin ciclotrón real (memoria permanente 3).",
        "Recomendación (no decisión operativa, 18.10): fijar strength=1.0 en corridas futuras y escalar "
        "explícitamente en el postproceso, para evitar ambigüedad de normalización.",
        "Valores de activación para el gemelo; no aptos para decisiones radiológicas reales.",
    ],
    "fuentes": ["evidencia_v072/fase70", "evidencia_v072/fase72", "f70_mc/full_voxel.json",
                "f70_mc/verif_voxel.json", "BON21", "ENDF/B-VIII.0 decay / NuDat",
                "Informe RadiopharmaOS v2.0 Digital Twins"],
    "md5_fuentes": {"full_voxel.json": md5(F70 / "full_voxel.json"),
                    "verif_voxel.json": md5(F70 / "verif_voxel.json"),
                    "diagnostico_ab": md5(OUT / "diagnostico_ab_strength.json")},
}
json.dump(out, open(OUT / "revision_unidades_tallies_f73.json", "w"), indent=1)

reg = OUT / "registro_unidades.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 73, "version": "2.9.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "revision_unidades_tallies_f73.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 73, "estado": "completada", "desv_rel_AB": desv_rel,
                  "A_sat_Bq_corregidas": {k: round(v["A_sat_Bq"], 1) for k, v in act.items()}}, indent=1))
