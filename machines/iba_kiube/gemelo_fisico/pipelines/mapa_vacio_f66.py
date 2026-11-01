#!/usr/bin/env python3
"""F66 — Escalón B: mapa de presión 3D por TPMC + pérdida de haz con P real.

Motor TPMC propio (vacuum_tpmc.py, algoritmo Molflow+ de TPMC-IPAC14/19 —
gitlab.cern.ch inaccesible desde el entorno; validación cruzada con el
binario oficial queda pendiente y declarada).

Cadena: cámara KIUBE (KIUBE-V04) → desgasificado paredes + H₂ de fuentes
(EST declarados) → 4 bombas ODP (s=0,25 EST) → mapa P(r) → gas_stripping
(F65) con perfil "mapa" → pérdida H⁻ con presión 3D, comparada con el
perfil uniforme/parametrizado del Escalón A.

Salida: evidencia_v072/fase66/mapa_vacio_f66.json + registro append-only.
"""
from __future__ import annotations

_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia externa
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

from app.domain.vacuum_tpmc import CamaraVacio, simular_mapa_presion, FUENTES_VACIO
from app.domain.gas_stripping import (PerfilPresion, simular_perdida_gas,
                                      FUENTES)

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase66")
EVID.mkdir(parents=True, exist_ok=True)

N_MOL = 600_000

# --- 1) mapa de presión TPMC ---
cam = CamaraVacio()  # cámara KIUBE-V04 por defecto
mapa = simular_mapa_presion(cam, n_moleculas=N_MOL, n_bins_r=20)

# --- 2) pérdida de haz con el mapa real vs perfiles del Escalón A ---
ESC = dict(corriente_ua=200.0, energia_final_mev=18.0,
           delta_e_kev_vuelta=98.5, b_media_t=1.360, r_pared_mm=548.0)

perfil_mapa = PerfilPresion.desde_mapa(mapa["r_mm"], mapa["p_mbar"])
p_media = sum(mapa["p_mbar"]) / len(mapa["p_mbar"])

casos = {
    "mapa_tpmc_f66": simular_perdida_gas(perfil=perfil_mapa, **ESC),
    "uniforme_media_mapa": simular_perdida_gas(
        perfil=PerfilPresion(p_base_mbar=p_media), **ESC),
    "uniforme_3e-7_referencia_f65": simular_perdida_gas(
        perfil=PerfilPresion(p_base_mbar=3e-7), **ESC),
}

resultados = {}
for nombre, r in casos.items():
    resultados[nombre] = {
        "fraccion_perdida": r.fraccion_perdida,
        "fraccion_perdida_pct": round(100 * r.fraccion_perdida, 6),
        "potencia_perdida_w": round(r.potencia_perdida_w, 4),
        "energia_media_h0_mev": round(r.parametros["energia_media_h0_mev"], 2),
        "impactos_azimut_36bins": [round(x, 6) for x in r.impactos_azimut],
    }

salida = {
    "fase": "F66",
    "nombre": "Mapa de presión TPMC (Escalón B) + pérdida H- con P 3D",
    "fecha_utc": datetime.now(timezone.utc).isoformat(),
    "estado": "simulacion_teorica_NO_produccion_real",
    "advertencia": ("Motor TPMC propio (algoritmo Molflow+ de literatura): "
                    "gitlab.cern.ch inaccesible (403) — validación cruzada con "
                    "Molflow+ oficial PENDIENTE. Flujos EST declarados: la "
                    "magnitud de P es tan buena como esos EST; la forma del "
                    "mapa es robusta. No cierra 18.19 ni PQ (14.9)."),
    "mapa_presion": mapa,
    "p_media_mbar": p_media,
    "escenario_haz": ESC,
    "resultados_perdida": resultados,
    "comparacion": {
        "mapa_vs_uniforme_misma_pmedia_pct": round(
            100 * (casos["mapa_tpmc_f66"].fraccion_perdida
                   / casos["uniforme_media_mapa"].fraccion_perdida - 1), 2),
    },
    "fuentes": {**FUENTES_VACIO, "stripping": FUENTES["secciones_stripping"]},
    "pendientes": [
        "validación cruzada con binario oficial Molflow+ (usuario sube zip)",
        "flujos de desgasificado/H₂ y s de bomba: medir o datasheets",
        "geometría 3D completa kiube (colinas 24 mm) vs gap efectivo",
        "acoplar mapa de impactos H⁰ al OpenMC de activación",
    ],
}

(EVID / "mapa_vacio_f66.json").write_text(
    json.dumps(salida, indent=2, ensure_ascii=False))

with (EVID / "registro_vacio.jsonl").open("a", encoding="utf-8") as f:
    f.write(json.dumps({
        "ts": salida["fecha_utc"], "fase": "F66", "estado": salida["estado"],
        "n_moleculas": N_MOL,
        "p_media_mbar": p_media,
        "perdida_mapa_pct": resultados["mapa_tpmc_f66"]["fraccion_perdida_pct"],
    }, ensure_ascii=False) + "\n")

print(f"P(r): centro {mapa['p_mbar'][0]:.2e} → pared {mapa['p_mbar'][-1]:.2e} mbar "
      f"(media {p_media:.2e}; Q/S = {mapa['p_equilibrio_q_s_mbar']:.2e})")
print(f"captura bombas: {mapa['absorbidos']}/{mapa['hits_bomba']} hits")
for k, v in resultados.items():
    print(f"{k}: pérdida {v['fraccion_perdida_pct']} %, {v['potencia_perdida_w']} W")
print("diferencia mapa vs uniforme:",
      salida["comparacion"]["mapa_vs_uniforme_misma_pmedia_pct"], "%")
print("OK F66 ->", EVID)
