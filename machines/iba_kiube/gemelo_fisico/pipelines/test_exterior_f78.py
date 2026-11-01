#!/usr/bin/env python3
"""Tests F78 — exterior F70 por extrapolación declarada (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, math

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase78")
E = json.load(open(EV / "exterior_f70_f78.json"))


def test_campos_canonicos():
    assert E["fase"] == 78 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert E["metodo"]["fuente_mu"].startswith("perfil de dosis F52")


def test_extrapolacion_declarada_y_fisica():
    m, r = E["metodo"], E["resultado"]
    assert m["corridas_directas"].startswith("13 corridas")
    assert "fallidas, declarado" in m["corridas_directas"]
    assert r["estado"].startswith("EXTRAPOLACION_DECLARADA")
    # μ positiva y en rango físico de concreto (0.005 - 0.05 cm^-1)
    assert 0.005 < m["mu_muro_cm^-1"] < 0.05
    # T estimada coherente con T directa F52 (mismo muro): factor <3
    t = r["T_H10_muro_este_estimada"] / r["T_referencia_F52_directa"]
    assert 1 / 3 < t < 3
    assert 0 < r["H10_losa_exterior_492_532_pSv_cm3_s_estimado"] < r["H10_losa_interior_282_292_pSv_cm3_s"]


def test_incertidumbre_propagada():
    assert E["resultado"]["factor_incertidumbre_1sigma"] > 1.0
    assert any("EXTRAPOLACIÓN DECLARADA" in a for a in E["advertencias"])
    assert any("pendiente_entorno" in a for a in E["advertencias"])


def test_dosis_sv_h_y_registro():
    d = E["dosis_exterior_Sv_h"]
    assert d["valor"] > 0 and d["estado"].startswith("estimado")
    reg = EV / "registro_exterior.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
