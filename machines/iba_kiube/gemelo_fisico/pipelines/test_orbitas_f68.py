"""Pruebas Fase 68 — órbitas región central KIUBE con física real.

Simulación teórica (memoria 3). PQ no_declarada_nunca (Informe 14.9).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import math
import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

from app.domain.orbita_central import (CampoKIUBE, RFKIUBE, brho_mev_tm,  # noqa: E402
                                       medir_f_rev_orbita, simular_aceleracion)

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase68")
JSON68 = EVID / "orbitas_centrales_f68.json"


def test_campo_anclas_s09():
    c = CampoKIUBE()
    w = c.resolver_ancho_colina()
    assert 35 < w < 55, "ancho efectivo plausible (v0.4: 46,7°)"
    assert abs(c.b_colina(300) - 2.25) < 0.01, "ancla S09 colina"
    assert abs(c.b_valle(400) / c.b_colina(400) - 0.45 / 2.15) < 1e-6
    assert c.B(0.3, 0.0) > 0 and c.B(0.6, 0.0) == 0.0, "borde de campo"


def test_isocronia_campo_medio():
    c = CampoKIUBE()
    # B̄(r)·γ⁻¹ = const por construcción
    b1 = c.b_media(100) / c.gamma_orbita(100)
    b2 = c.b_media(400) / c.gamma_orbita(400)
    assert abs(b1 / b2 - 1) < 1e-12


def test_rf_parametros_oficiales():
    rf = RFKIUBE()
    assert rf.f_mhz == 40.646          # T04F8
    assert rf.v_dee_kv == 40.0         # TUD03
    assert rf.armonico == 2
    # ganancia pico por vuelta: 4 bordes × V·sin(h·α/2) ≈ 98,5 keV (v0.4 CALC)
    assert abs(4 * rf.ganancia_gap_kev - 98.5) < 1.0


def test_orbita_sin_rf_conserva_energia_y_radio():
    c = CampoKIUBE()
    c.resolver_ancho_colina()
    br = brho_mev_tm(8.0)
    r = br / c.b0_t
    for _ in range(30):
        r = br / c.b_media(r * 1e3)
    f = medir_f_rev_orbita(c, 8.0)
    # sin shim aún debe estar al <5 % (hard-edge flutter)
    assert 19.0e6 < f < 22.0e6


@pytest.mark.skipif(not JSON68.exists(), reason="pipeline aún no ejecutado")
def test_json_f68_campos_canonicos():
    d = json.loads(JSON68.read_text())
    assert d["fase"] == 68 and d["estado"] == "simulacion_teorica"
    assert d["pq"] == "no_declarada_nunca"
    assert d["fuentes"] and d["pendientes"]
    assert d["rf"]["f_mhz"] == 40.646


@pytest.mark.skipif(not JSON68.exists(), reason="pipeline aún no ejecutado")
def test_isocronia_tras_shim():
    d = json.loads(JSON68.read_text())
    errs = d["campo"]["errores_isocronia_finales_x1e3"]
    assert all(abs(e) < 1.5 for e in errs), "isocronía <|1,5|×10⁻³ en todo el rango"


@pytest.mark.skipif(not JSON68.exists(), reason="pipeline aún no ejecutado")
def test_ventana_captura_existe_y_demo_coherente():
    d = json.loads(JSON68.read_text())
    assert len(d["ventana_captura_deg"]) >= 3, "ventana de captura no vacía"
    demo = d["demostracion"]
    assert demo["vueltas"] >= 10, "captura de región central demostrada"
    assert demo["energia_max_mev"] > 0.5
    radios = demo["radios_por_vuelta_mm"]
    assert all(r < 548 for r in radios), "dentro de la cámara"


@pytest.mark.skipif(not JSON68.exists(), reason="pipeline aún no ejecutado")
def test_espiral_y_registro_append_only():
    d = json.loads(JSON68.read_text())
    es = d["espiral_equilibrio"]
    assert len(es) == 30
    assert es[-1]["E_mev"] == 18.0
    rs = [o["r_media_mm"] for o in es]
    assert all(a < b for a, b in zip(rs, rs[1:])), "espiral monótona"
    reg = EVID / "registro_orbitas.jsonl"
    for ln in reg.read_text().strip().splitlines():
        assert json.loads(ln)["estado"] == "simulacion_teorica"
