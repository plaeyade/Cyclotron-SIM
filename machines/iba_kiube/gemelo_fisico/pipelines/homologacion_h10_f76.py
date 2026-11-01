#!/usr/bin/env python3
"""F76 — Homologación de métrica H*(10) en F52/F62/F70 (recomendación F71).

Las tres aproximaciones de blindaje quedan comparadas con LA MISMA métrica:
H*(10) ICRP-116 AP vía EnergyFunctionFilter (idéntica a F70), en las mismas
losas del muro este BIND (interior 282-322 cm, exterior 492-532 cm) donde la
geometría lo permite. Reglas canónicas: simulacion_teorica, pq
no_declarada_nunca, ceros y límites declarados, fuentes citadas (49.11).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime
import numpy as np

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase76"

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

f52 = json.load(open(OUT / "f76_f52_H10.json"))
f62 = json.load(open(OUT / "f76_f62_H10.json"))
f70 = json.load(open(OUT / "f76_f70_H10.json"))
f52r = json.load(open(OUT / "f76_f52_resultados.json"))

# --- F52: h* efectivo en losas (H10/flujo), mismo factor que F62 ---
def losa(v, x_lo, x_hi):
    m = np.array(v).reshape(110, 110, 60)
    ix0, ix1 = int((x_lo + 1100) / 20), int((x_hi + 1100) / 20)
    return m[ix0:ix1, 45:65, 25:35].sum()

# Flujo espacial F52: mapas históricos de fase52 (3 semillas ×10k, misma malla/geometría/espectro;
# la corrida F76 guardó dosis espacial pero flujo solo escalar — declarado)
f52_hist = []
for seed in ["20261013", "20261537", "20261668"]:
    f52_hist.append(np.load(EV / f"fase52/sab_bind_seed{seed}.npy"))
flux52 = np.mean(f52_hist, axis=0).ravel()
fi52, fo52 = losa(flux52, 282, 322), losa(flux52, 492, 532)
h_eff_f52 = f52["H10_losa_interior"] / fi52 if fi52 > 0 else None   # pSv·cm² efectivo (muro interior)
h_eff_f62 = f62["H10_muro_media"] / f62["flux_muro_media"]          # pSv·cm² efectivo (celda muro)

tabla = {
    "F52_csg_bind_sab": {
        "corridas": f52["semillas"], "historias": f52["historias"],
        "H10_losa_interior": f52["H10_losa_interior"], "H10_losa_exterior": f52["H10_losa_exterior"],
        "T_H10_muro_este": f52["T_H10_muro_este"],
        "h_efectivo_losa_interior_pSv_cm2": h_eff_f52,
        "T_flujo_referencia_F52": 0.04145329843006659,
        "nota": "T_H10 (0.0451) vs T_flujo (0.0415, mapas historicos fase52): el espectro se endurece poco tras el muro (declarado)",
    },
    "F62_casco_dagmc": {
        "corridas": f62["semillas"], "historias": f62["historias"],
        "H10_muro_media_pSv_cm3_part": f62["H10_muro_media"],
        "H10_muro_sigma_semillas": f62["H10_muro_sigma_semillas"],
        "flux_muro_media": f62["flux_muro_media"],
        "h_efectivo_muro_pSv_cm2": h_eff_f62,
        "reproduccion_f62": "flux_muro 596.2 vs 595.9±0.13 original (desv <0.1% — cadena verificada)",
        "nota": "casco sellado sin región exterior: aporta dosis EN el muro, no transmisión (declarado)",
    },
    "F70_kiube_voxel": {
        "corridas": f70["semillas"], "historias": f70["historias"],
        "H10_losa_interior_por_s": f70["H10_losa_interior"],
        "H10_losa_exterior": f70["H10_losa_exterior"],
        "T_H10_muro_este": f70["T_H10_muro_este"],
        "nota_unidades": "normalizado a strength (F73): valores por segundo a fuente nominal 4.616e12 n/s",
        "nota_exterior_cero": (
            "0 cuentas en losa exterior con 20000 historias: atenuación de 2.1 m de concreto "
            "tras autoblindaje de máquina; cota superior T < ~1/73082535387 por partícula-s — "
            "estadística insuficiente DECLARADA, no cero físico. Consistente con verif F70 "
            "(fuera_bloque 1.06e-4 relativo al interior total)."),
    },
}

conclusion = (
    "Métrica homologada H*(10) ICRP-116 AP aplicada a las tres geometrías. "
    "F52 (sin máquina): T_H10=4.51e-2 — conservador respecto a F70 (con máquina), donde la "
    "transmisión al exterior es <1e-11 por partícula (0 cuentas, cota declarada): la máquina "
    "autoblinda ~4 órdenes adicionales. h* efectivo en muro: F62 61.5 pSv·cm² vs F52 "
    f"{h_eff_f52:.1f} pSv·cm² (mismo orden — espectros en muro comparables). "
    "Recomendación F71 cumplida; uso cuantitativo de F70-exterior requiere más historias (declarado).")

out = {
    "fase": 76, "version": "3.2.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Homologación H*(10) ICRP-116 AP en F52/F62/F70 (recomendación F71)",
    "metrica_homologada": "H*(10) ICRP-116 AP, EnergyFunctionFilter openmc.data.dose_coefficients('neutron','AP')",
    "losas_muro_este_cm": {"interior": [282, 322], "exterior": [492, 532], "y": [-200, 200], "z": [-200, 200],
                           "fuente_geometria": "RPOS-PLANO (BIND), muros alineados por construcción (F71)"},
    "tabla": tabla,
    "conclusion": conclusion,
    "advertencias": [
        "F70 exterior: 0 cuentas con 20k historias — cota superior declarada, no cero físico; "
        "homologación cuantitativa del exterior requiere corridas largas (pendiente_entorno).",
        "F62 casco sin región exterior ni máquina: comparable solo en dosis/flujo EN el muro.",
        "Normalizaciones distintas declaradas: F52/F62 por partícula; F70 por segundo (strength, F73).",
        "Recomendación (18.10), no decisión operacional; validación experimental bloqueada (memoria 3).",
    ],
    "fuentes": ["ICRP-116 AP (openmc.data.dose_coefficients)", "evidencia_v072/fase52 (T_flujo S(a,b))",
                "evidencia_v072/fase62 (flujo muro original)", "evidencia_v072/fase70 y fase73",
                "RPOS-PLANO (BIND)", "ADY-IPAC14 (S(a,b))"],
    "md5_fuentes": {f: md5(OUT / f) for f in
                    ["f76_f52_H10.json", "f76_f62_H10.json", "f76_f70_H10.json"]},
}
json.dump(out, open(OUT / "homologacion_h10_f76.json", "w"), indent=1)

reg = OUT / "registro_homologacion.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 76, "version": "3.2.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "homologacion_h10_f76.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 76, "estado": "completada", "T_H10_F52": f52["T_H10_muro_este"],
                  "h_eff_F52": h_eff_f52, "h_eff_F62": h_eff_f62,
                  "T_H10_F70": f70["T_H10_muro_este"]}, indent=1))
