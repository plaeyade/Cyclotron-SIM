#!/usr/bin/env python3
"""Tests F65 — pérdida por stripping en gas residual."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia externa
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

from app.domain.gas_stripping import (PerfilPresion, beta2, radio_orbita_mm,
                                      sigma_mezcla_cm2, simular_perdida_gas,
                                      COMPOSICION_RESIDUAL_NAKHAI)

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase65")


def test_sigma_decrece_con_energia():
    # física: σ = C/β² decrece con E (PLUM-CYR)
    s_lo = sigma_mezcla_cm2(COMPOSICION_RESIDUAL_NAKHAI, 1.0)
    s_hi = sigma_mezcla_cm2(COMPOSICION_RESIDUAL_NAKHAI, 18.0)
    assert s_lo > s_hi > 0


def test_radio_orbita_fisico():
    # r(18 MeV, 1.36 T) ≈ 453 mm < 548 mm de cámara (KIUBE-V04)
    assert 440 < radio_orbita_mm(18.0, 1.36) < 470
    assert radio_orbita_mm(0.1, 1.36) < radio_orbita_mm(18.0, 1.36)


def test_perdida_monotona_en_presion():
    kw = dict(corriente_ua=200, energia_final_mev=18.0,
              delta_e_kev_vuelta=98.5, b_media_t=1.360)
    lo = simular_perdida_gas(perfil=PerfilPresion(1e-7), **kw)
    hi = simular_perdida_gas(perfil=PerfilPresion(1e-6), **kw)
    assert 0 < lo.fraccion_perdida < hi.fraccion_perdida < 1
    # proporcionalidad casi lineal en régimen de pérdida pequeña
    assert abs(hi.fraccion_perdida / lo.fraccion_perdida - 10.0) < 0.5


def test_perdida_en_rango_fisico():
    # nominal 3e-7 mbar: pérdida del orden de 1e-3 – 1e-1 (0,1–10 %). F82: hasta v3.7.0 el rango esperado era
    # 1e-5 – 1e-3 porque MBAR_A_PA valía 0,1 (1 mbar = 100 Pa); la evidencia histórica de F65 no se reescribe.
    r = simular_perdida_gas(corriente_ua=200, energia_final_mev=18.0,
                            delta_e_kev_vuelta=98.5, b_media_t=1.360,
                            perfil=PerfilPresion(3e-7))
    assert 1e-3 < r.fraccion_perdida < 1e-1
    assert r.potencia_perdida_w > 0


def test_impactos_normalizados_y_en_pared():
    r = simular_perdida_gas(corriente_ua=200, energia_final_mev=18.0,
                            delta_e_kev_vuelta=98.5, b_media_t=1.360,
                            perfil=PerfilPresion(3e-7))
    assert len(r.impactos_azimut) == 36
    assert abs(sum(r.impactos_azimut) - 1.0) < 1e-9
    assert all(x >= 0 for x in r.impactos_azimut)


def test_perfil_radial_mayor_que_base_en_centro():
    p = PerfilPresion(p_base_mbar=3e-7, modelo="radial", p_centro_mbar=1e-6)
    assert p.presion_mbar(0.0) == 1e-6
    assert p.presion_mbar(548.0) == 3e-7
    assert 3e-7 < p.presion_mbar(200.0) < 1e-6


def test_f65_json_coherente():
    d = json.loads((EVID / "perdida_gas_f65.json").read_text())
    assert "NO_produccion_real" in d["estado"]
    assert len(d["resultados"]) == 4
    assert "PLUM-CYR" in d["fuentes"]["secciones_stripping"]
    assert "Molflow" in " ".join(d["pendientes"])
    nom = d["resultados"]["vacio_nominal_3e-7"]
    assert 1e-5 < nom["fraccion_perdida"] < 1e-3


def test_registro_append_only():
    reg = EVID / "registro_perdida_gas.jsonl"
    lineas = [l for l in reg.read_text().strip().splitlines() if l]
    assert len(lineas) >= 1
    for l in lineas:
        e = json.loads(l)
        assert e["fase"] == "F65" and "casos" in e
