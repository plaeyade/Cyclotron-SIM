"""Tests Fase 69 — Recalibrado del campo con mapa medido KEY22.

Verifica: carga del CSV medido (huecos interpolados, declarados), anclas del
perfil (⟨B⟩(450 mm) y joroba central vs publicado), descomposición por flutter
medido, convergencia del shim radial, JSON canónico (49.11, 14.9, 18.10),
órbita 18 MeV coherente con F68, registro append-only y compatibilidad F68
(sin KEY22 el campo es el analítico γ·B0 de siempre).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import math
import sys
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)

from app.domain.orbita_central import (  # noqa: E402
    CampoKIUBE, RFKIUBE, cargar_perfil_key22, calibrar_shim_radial,
    medir_f_rev_orbita, brho_mev_tm)

CSV = RAIZ / "data" / "kiube" / "kiube_Bavg_flutter_KEY22.csv"
EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase69")


def test_csv_existe_y_carga():
    assert CSV.exists()
    c = cargar_perfil_key22(CampoKIUBE(), CSV)
    assert len(c.key22_r_mm) == 506
    assert min(c.key22_r_mm) == pytest.approx(5.8, abs=0.2)
    assert max(c.key22_r_mm) == pytest.approx(549.1, abs=0.2)
    assert c.gap_valle_mm == 110.0  # KEY22-GEOM
    # huecos interpolados: sin NaN en ningún punto
    assert not any(math.isnan(b) for b in c.key22_bavg_t)
    assert not any(math.isnan(f) for f in c.key22_f)


def test_perfil_medido_anclas():
    c = cargar_perfil_key22(CampoKIUBE(), CSV)
    # ⟨B⟩(450 mm) = 1,3316 T (KEY22 Bavg_modelo 13315,5 G; medido IBA 13252,7 G)
    assert c.b_media(450.0) == pytest.approx(1.3316, abs=0.02)
    # joroba central medida ≈ 1,033 (KEY22 B_centro_rel = 1,0329)
    joroba = c.b_media(23.0) / c.b_media(450.0)
    assert joroba == pytest.approx(1.033, abs=0.01)


def test_flutter_medido_colina_valle():
    c = cargar_perfil_key22(CampoKIUBE(), CSV)
    bm, bh, bv = c.b_media(450.0), c.b_colina(450.0), c.b_valle(450.0)
    assert bh > bm > bv
    # coherente con FE kiube v0.7 (2,169/0,482 T; definición distinta, ±8 %)
    assert bh == pytest.approx(2.17, abs=0.12)
    assert bv == pytest.approx(0.50, abs=0.10)


def test_shim_radial_converge_sobre_medida():
    c = cargar_perfil_key22(CampoKIUBE(), CSV)
    rf = RFKIUBE()
    info = calibrar_shim_radial(c, rf, energias_mev=(0.1, 1.0, 8.0, 18.0),
                                iteraciones=3)
    assert all(abs(e) < 3e-3 for e in info["errores_finales"])
    f_obj = rf.f_mhz * 1e6 / rf.armonico
    assert medir_f_rev_orbita(c, 8.0) == pytest.approx(f_obj, rel=3e-3)


def test_json_canonico():
    d = json.loads((EVID / "orbitas_recal_f69.json").read_text())
    assert d["fase"] == 69 and d["version"] == "2.5.0"
    assert d["estado"] == "simulacion_teorica"
    assert d["pq"] == "no_declarada_nunca"                    # 14.9
    assert d["campo"]["fuente_perfil"] == "KEY22"             # 49.11
    assert "KEY22" in d["fuentes"] and "S09" in d["fuentes"]
    assert d["campo"]["error_isocronia_denso_60ptos_max"] < 1.5e-3
    assert d["captura_rf_mejor_fase"]["vueltas"] >= 3
    assert "pendiente" in d["captura_rf_mejor_fase"]["nota"]


def test_orbita_18mev_coherente_con_f68():
    d = json.loads((EVID / "orbitas_recal_f69.json").read_text())
    o69, o68 = d["orbita_18mev"]["f69_key22"], d["orbita_18mev"]["f68_analitico"]
    for o in (o69, o68):
        assert 380 < o["r_min_mm"] < o["r_max_mm"] < 512
    assert abs(o69["r_min_mm"] - o68["r_min_mm"]) / o68["r_min_mm"] < 0.05
    assert abs(o69["r_max_mm"] - o68["r_max_mm"]) / o68["r_max_mm"] < 0.05


def test_registro_append_only():
    reg = EVID / "registro_orbitas.jsonl"
    lineas = [json.loads(l) for l in reg.read_text().splitlines() if l.strip()]
    fases = [l.get("fase") for l in lineas]
    assert 69 in fases
    assert len(lineas) == len(set(l["ts_utc"] for l in lineas))  # sin reescribir


def test_compatibilidad_f68_sin_key22():
    c = CampoKIUBE()
    assert not c.key22_bavg_t
    # sin perfil KEY22: b_media es el analítico γ·B0 de F68
    esperado = c.gamma_orbita(300.0) * c.b0_t
    assert c.b_media(300.0) == pytest.approx(esperado, rel=1e-12)
    assert c.b_valle(300.0) == pytest.approx(
        c.rho_vh * c.b_colina(300.0), rel=1e-12)
