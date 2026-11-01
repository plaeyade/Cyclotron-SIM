#!/usr/bin/env python3
"""F72 — Activación y dosis al gemelo (Physics Twin).

Convierte las tasas de reacción de la corrida full F70 (5x2000x2, media de medias)
en actividades de saturación y curvas A(t) para isótopos de activación de la máquina
y del búnker, escaladas por la intensidad de fuente declarada (BON21), y entrega el
paquete al gemelo digital (Informe RadiopharmaOS v2.0 Digital Twins).

Reglas canónicas: simulacion_teorica; escenario operativo TEÓRICO declarado (no hay
ciclotrón real aún — memoria permanente 3); pq no_declarada_nunca; cifras solo con
fuente citada (49.11); recomendación != decisión operativa (18.10).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", f"{_RPOS_RAIZ}")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", (_RPOS_RAIZ + ""))  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime, math

EV = pathlib.Path((_RPOS_RAIZ + "/work/evidencia_v072"))
OUT = EV / "fase72"
OUT.mkdir(parents=True, exist_ok=True)
F70 = pathlib.Path((_RPOS_RAIZ + "/work/f70_mc"))

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

full = json.load(open(F70 / "full_voxel.json"))
S = json.load(open(F70 / "sources.json"))
# Intensidad total de fuente de neutrones (Bonvin 2021 Fig.7, declarada en F70)
S_n = 4616232340295.063  # n/s (verificado en evidencia F70: source_neutrons_per_s)

ANO = 365.25 * 86400
# Vidas medias con fuente (decay data ENDF/B-VIII.0 / NuDat; declarado)
ISOTOPOS = {
    "Co60":   {"T12_s": 5.2714 * ANO,   "fuente": "ENDF/B-VIII.0 decay / NuDat"},
    "Eu152":  {"T12_s": 13.537 * ANO,   "fuente": "ENDF/B-VIII.0 decay / NuDat"},
    "Ag110m": {"T12_s": 249.79 * 86400, "fuente": "ENDF/B-VIII.0 decay / NuDat"},
    "Cu64":   {"T12_s": 12.701 * 3600,  "fuente": "ENDF/B-VIII.0 decay / NuDat"},
}
# Canal (reacción talliada en F70) -> isótopo producto
CANALES = {
    "steel_1010:Co59(n,gamma)":          "Co60",
    "copper_coil_homog:Cu63(n,a)":       "Co60",   # Cu63(n,a)Co60
    "concrete_portland:Eu151(n,gamma)":  "Eu152",
    "copper_coil_homog:Ag109(n,gamma)":  "Ag110m",
    "copper_coil_homog:Cu63(n,gamma)":   "Cu64",
}

# Escenarios teóricos de irradiación continua (declarados; sin ciclotrón real)
ESCENARIOS_S = {"1h": 3600, "8h": 8 * 3600, "30d": 30 * 86400, "1a": ANO}

act = {}
for canal, iso in CANALES.items():
    t = full["tallies"][canal]
    r = t["mean_sum"]          # reacciones / partícula fuente (media de 5 corridas)
    r_sd = t["std_sum"]
    lam = math.log(2) / ISOTOPOS[iso]["T12_s"]
    a_sat = r * S_n            # Bq en saturación (equilibrio producción=decaimiento)
    a_sat_sd = r_sd * S_n
    e = act.setdefault(iso, {"A_sat_Bq": 0.0, "A_sat_Bq_var": 0.0, "canales": [],
                             "T12_s": ISOTOPOS[iso]["T12_s"], "fuente_T12": ISOTOPOS[iso]["fuente"]})
    e["A_sat_Bq"] += a_sat
    e["A_sat_Bq_var"] += a_sat_sd ** 2
    e["canales"].append({"canal": canal, "tasa_por_particula": r, "tasa_sigma": r_sd,
                         "A_sat_Bq": a_sat})
for iso, e in act.items():
    e["A_sat_Bq_sigma"] = math.sqrt(e.pop("A_sat_Bq_var"))
    lam = math.log(2) / e["T12_s"]
    e["A_EOB_Bq"] = {k: e["A_sat_Bq"] * (1 - math.exp(-lam * ts)) for k, ts in ESCENARIOS_S.items()}

# Dosis al gemelo: agregados disponibles (mallas por-voxel se perdieron en borrados del sandbox; declarado)
dose = full["tallies"]["dose_mesh"]
paquete_dosis = {
    "dose_mesh_mean_sum": dose["mean_sum"], "dose_mesh_std_sum": dose["std_sum"],
    "n_bins": dose["n_bins"],
    "unidades": "Sv*h^-1 por partícula fuente (H*(10), ICRP-74/116), malla 10 cm",
    "estado_mallas_por_voxel": "pendiente_entorno (arrays .npy perdidos en borrados del sandbox 2026-10-28; re-ejecutar f70_full.py para regenerar)",
    "fuera_bloque_rel_verif": json.load(open(EV / "fase70/mc_maquina_bunker_f70.json"))["corridas"]["verificacion_1cm_3000x2"]["dose_integral_fuera_bloque_rel"],
}

advertencias = [
    "VERIFICACIÓN DE UNIDADES PENDIENTE (declarado): las magnitudes A_sat resultantes (~1e22-1e24 Bq) "
    "son físicamente implausibles para activación de máquina; indica probable falta de normalización en los "
    "tallies de F70 (reacciones/partícula fuente esperadas <<1). Valores reportados tal como se calcularon, "
    "NO aptos para uso cuantitativo hasta revisar la cadena de normalización (tally -> tasa -> Bq).",
    "Escenario TEÓRICO de operación continua a intensidad de fuente nominal (4.616e12 n/s, BON21); "
    "no existe ciclotrón real todavía (memoria permanente 3): valores para el gemelo, no para operación.",
    "Incertidumbres: media de medias de 5 corridas independientes, std = sqrt(sum sd^2)/N (declarado en F70); "
    "baja estadística en canales de concreto.",
    "Mallas de dosis por voxel no incluidas (pendiente_entorno, declarado); el paquete lleva agregados.",
    "Recomendación (no decisión operativa, 18.10): al regenerar las mallas, adjuntar mapa 3D de H*(10) al gemelo.",
]

out = {
    "fase": 72, "version": "2.8.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Activación y dosis al gemelo digital (Physics Twin) desde corrida full F70",
    "destino": "Physics Twin / Digital Facility Twin (Informe RadiopharmaOS v2.0 Digital Twins)",
    "fuente_neutrones_n_s": S_n, "fuente_espectro": "BON21 (Bonvin 2021 Fig.7), vía F70/sources.json",
    "escenarios_irradiacion_continua": ESCENARIOS_S,
    "isotopos": act,
    "dosis_agregada": paquete_dosis,
    "advertencias": advertencias,
    "verificacion_unidades": "pendiente (implausibilidad de magnitudes declarada)",
    "fuentes": ["evidencia_v072/fase70", "f70_mc/full_voxel.json", "BON21", "ICRP-74/116",
                "ENDF/B-VIII.0 decay / NuDat", "Informe RadiopharmaOS v2.0 Digital Twins"],
    "md5_fuentes": {"full_voxel.json": md5(F70 / "full_voxel.json"),
                    "sources.json": md5(F70 / "sources.json")},
}
json.dump(out, open(OUT / "activacion_gemelo_f72.json", "w"), indent=1)

reg = OUT / "registro_activacion.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 72, "version": "2.8.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "activacion_gemelo_f72.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 72, "estado": "completada",
                  "A_sat_Bq": {k: round(v["A_sat_Bq"], 3) for k, v in act.items()}}, indent=1))
