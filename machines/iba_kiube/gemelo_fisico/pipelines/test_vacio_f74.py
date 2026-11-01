#!/usr/bin/env python3
"""Tests F74 — validación Chiggiato del vacío (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase74")
E = json.load(open(EV / "vacio_chiggiato_f74.json"))


def test_campos_canonicos():
    assert E["fase"] == 74 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert any("CHIGGIATO-CAS" in f for f in E["fuentes"])
    assert any("CLAUSING" in f for f in E["fuentes"])


def test_consistencia_triple_s_efectiva():
    s = E["comparacion_triple"]["S_eff_ls"]
    # apertura ideal Chiggiato vs TPMC y Molflow: acuerdo <1%
    assert abs(s["chiggiato_apertura_ideal"] - s["tpmc_f66"]) / s["tpmc_f66"] < 0.01
    assert abs(s["chiggiato_apertura_ideal"] - s["molflow_f67"]) / s["molflow_f67"] < 0.01
    # corrección de canal corto reduce S (física) pero no más de 50%
    assert 0.5 < s["chiggiato_canal_corto"] / s["chiggiato_apertura_ideal"] < 1.0


def test_p_eq_coherente():
    p = E["comparacion_triple"]["P_eq_mbar"]
    # P=Q/S ideal coincide con TPMC (<1%)
    assert abs(p["chiggiato_ideal"] - p["tpmc_f66"]) / p["tpmc_f66"] < 0.01
    # canal corto da P mayor que ideal
    assert p["chiggiato_canal_corto"] > p["chiggiato_ideal"]
    # dentro del rango físico del mapa TPMC (orden 1e-7)
    assert 5e-8 < p["chiggiato_canal_corto"] < 1e-6


def test_honestidad_y_registro():
    assert any("Clausing" in a for a in E["advertencias"])
    assert any("BLOQUEADA hasta haz real" in a for a in E["advertencias"])
    assert E["comparacion_triple"]["W_clausing_puerto"] == round(1 / (1 + 3 * 132 / (8 * 127)), 4)
    reg = EV / "registro_vacio.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
