#!/usr/bin/env python3
"""F65 — Pérdida de haz H⁻ por stripping en gas residual (Escalón A).

Escenario: Cyclone KIUBE 180 (doble haz, 2 × 100 µA — KIUBE-OM Tabla 1-1),
órbita del modelo kiube v0.4 (ΔE = 98,5 keV/vuelta, <B> = 1,360 T,
cámara r = 548 mm — KIUBE-V04), gas residual 70 % H₂O / 30 % aire
(NAKHAI2022), secciones σ = C/β² (PLUM-CYR / Gillespie).

Perfil de presión PARAMETRIZADO y DECLARADO (supuesto_declarado_pendiente_
fuente): rango 1×10⁻⁸ – 1×10⁻⁶ mbar (vacío típico de ciclotrón compacto
declarado en literatura; el valor de diseño del KIUBE no es público).
El mapa 3D real lo aportará Molflow+ en el Escalón B.

Salida: evidencia_v072/fase65/perdida_gas_f65.json + registro append-only.
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

from app.domain.gas_stripping import (FUENTES, PerfilPresion,
                                      COMPOSICION_RESIDUAL_NAKHAI,
                                      simular_perdida_gas)

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase65")
EVID.mkdir(parents=True, exist_ok=True)

# Escenario KIUBE 180 dual beam (KIUBE-OM Tabla 1-1; kiube v0.4 params.yaml)
ESCENARIO = {
    "maquina": "Cyclone KIUBE 180 (doble haz 2 × 100 µA)",
    "corriente_ua": 200.0,            # total en 2 blancos
    "energia_final_mev": 18.0,        # KIUBE-V04 (TUP03); 17,5 en blanco (KIUBE-OM)
    "delta_e_kev_vuelta": 98.5,       # KIUBE-V04 (CALC)
    "b_media_t": 1.360,               # KIUBE-V04 (CALC)
    "r_pared_mm": 548.0,              # KIUBE-V04 (T04F2, cámara)
}

# Barrido de presión declarado: base y pico central (fuentes de iones H₂)
CASOS = {
    "vacio_buenp_1e-7":  PerfilPresion(p_base_mbar=1e-7, modelo="uniforme"),
    "vacio_nominal_3e-7": PerfilPresion(p_base_mbar=3e-7, modelo="uniforme"),
    "vacio_degradado_1e-6": PerfilPresion(p_base_mbar=1e-6, modelo="uniforme"),
    "perfil_radial_3e-7_1e-6": PerfilPresion(
        p_base_mbar=3e-7, modelo="radial", p_centro_mbar=1e-6),
}

resultados = {}
for nombre, perfil in CASOS.items():
    r = simular_perdida_gas(
        corriente_ua=ESCENARIO["corriente_ua"],
        energia_final_mev=ESCENARIO["energia_final_mev"],
        delta_e_kev_vuelta=ESCENARIO["delta_e_kev_vuelta"],
        b_media_t=ESCENARIO["b_media_t"],
        perfil=perfil,
        r_pared_mm=ESCENARIO["r_pared_mm"],
    )
    resultados[nombre] = {
        "fraccion_perdida": r.fraccion_perdida,
        "fraccion_perdida_pct": round(100 * r.fraccion_perdida, 4),
        "potencia_perdida_w": round(r.potencia_perdida_w, 4),
        "n_vueltas": r.n_vueltas,
        "energia_media_h0_mev": round(r.parametros["energia_media_h0_mev"], 2),
        "impactos_azimut_36bins": [round(x, 6) for x in r.impactos_azimut],
        "perdida_vuelta_max": max(r.perdida_por_vuelta),
        "perdida_vuelta_min": min(r.perdida_por_vuelta),
        "vuelta_de_max_perdida": r.perdida_por_vuelta.index(max(r.perdida_por_vuelta)) + 1,
    }

salida = {
    "fase": "F65",
    "nombre": "Pérdida de haz H- por stripping en gas residual (Escalón A)",
    "fecha_utc": datetime.now(timezone.utc).isoformat(),
    "estado": "simulacion_teorica_NO_produccion_real",
    "advertencia": ("Modelo teórico con P(r) parametrizada declarada; no usa "
                    "datos de haz reales (memoria 3). No cierra 18.19 ni PQ (14.9)."),
    "escenario": ESCENARIO,
    "composicion_gas": COMPOSICION_RESIDUAL_NAKHAI,
    "resultados": resultados,
    "fuentes": FUENTES,
    "pendientes": [
        "P(r) parametrizado → mapa Molflow+ 3D (Escalón B)",
        "σ(E) exactas por gas → tablas Gillespie completas",
        "impactos 2D plano medio → trazado 3D con geometría kiube",
        "alimentar término fuente distribuido al OpenMC de activación (F66+)",
    ],
}

(EVID / "perdida_gas_f65.json").write_text(
    json.dumps(salida, indent=2, ensure_ascii=False))

registro = EVID / "registro_perdida_gas.jsonl"
with registro.open("a", encoding="utf-8") as f:  # append-only ALCOA+
    f.write(json.dumps({
        "ts": salida["fecha_utc"], "fase": "F65", "estado": salida["estado"],
        "casos": {k: v["fraccion_perdida_pct"] for k, v in resultados.items()},
    }, ensure_ascii=False) + "\n")

for k, v in resultados.items():
    print(f"{k}: pérdida {v['fraccion_perdida_pct']} %, "
          f"{v['potencia_perdida_w']} W, E_H0 ~{v['energia_media_h0_mev']} MeV, "
          f"máx pérdida en vuelta {v['vuelta_de_max_perdida']}")
print("OK F65 ->", EVID)
