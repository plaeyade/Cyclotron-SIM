#!/usr/bin/env python3
"""Tests F73 — revisión de unidades de tallies F70 (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase73")
E = json.load(open(EV / "revision_unidades_tallies_f73.json"))
S_n = 4616232340295.063


def test_campos_canonicos():
    assert E["fase"] == 73 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"


def test_prueba_ab_consistente():
    h = E["hallazgo"]
    # la razón medida debe coincidir con la fuerza de fuente (<1% desviación)
    assert h["desviacion_relativa"] < 0.01
    assert abs(h["razon_esperada_strength"] - 2308116170147.5317) < 1e6


def test_activacion_corregida_plausible():
    # corregidas deben ser las de F72 divididas por S_n
    f72 = json.load(open(f"{_RPOS_RAIZ}/work/evidencia_v072/fase72/activacion_gemelo_f72.json"))
    for iso, v in E["activacion_corregida"].items():
        esperado = f72["isotopos"][iso]["A_sat_Bq"] / S_n
        assert abs(v["A_sat_Bq"] - esperado) / esperado < 1e-9
        assert v["A_sat_Bq"] > 0
        a = v["A_EOB_Bq"]
        assert a["1h"] <= a["8h"] <= a["30d"] <= a["1a"] <= v["A_sat_Bq"] * 1.0000001
    # magnitudes ya no absurdas para activación de máquina (< 1e15 Bq)
    assert all(v["A_sat_Bq"] < 1e15 for v in E["activacion_corregida"].values())


def test_correcciones_historicas_y_registro():
    ch = E["correcciones_historicas"]
    assert ch["verif_max_dose_rate_Sv_h"]["corregido"] < ch["verif_max_dose_rate_Sv_h"]["original"]
    assert ch["verif_max_dose_rate_Sv_h"]["corregido"] > 0
    assert E["dosis_agregada_corregida"]["fuera_bloque_rel_verif"] < 0.05
    assert "append-only" in E["politica"]
    reg = EV / "registro_unidades.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
