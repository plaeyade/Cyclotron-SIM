#!/usr/bin/env python3
"""Tests F79 — propiedades termo-mecánicas con fuente (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase79")
E = json.load(open(EV / "propiedades_fem_f79.json"))
MATS = json.load(open(f"{_RPOS_RAIZ}/work/f70_mc/materials.json"))


def test_campos_canonicos():
    assert E["fase"] == 79 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert len(E["propiedades"]) == len(MATS)


def test_todo_valor_tiene_fuente_o_estado_honesto():
    for nombre, p in E["propiedades"].items():
        est = p["estado"]
        assert est.startswith(("con_fuente", "calc_declarado", "pendiente_fuente", "no_aplica"))
        if est == "con_fuente":
            assert p["fuente"] and all(p[k] is not None for k in ("k_W_mK", "cp_J_kgK", "E_GPa", "nu"))
        if est == "pendiente_fuente":
            assert all(p[k] is None for k in ("k_W_mK", "cp_J_kgK", "E_GPa", "nu"))


def test_valores_en_rango_fisico():
    for p in E["propiedades"].values():
        if p["estado"] in ("con_fuente", "calc_declarado"):
            assert 0.05 < p["k_W_mK"] < 500
            assert 100 < p["cp_J_kgK"] < 2000
            assert 1 < p["E_GPa"] < 300
            assert 0.1 < p["nu"] < 0.5
    # anclas conocidas
    assert E["propiedades"]["steel_1010"]["k_W_mK"] == 49.8
    assert E["propiedades"]["stainless_304L"]["E_GPa"] == 193
    assert E["propiedades"]["copper_ofe"]["nu"] == 0.34


def test_bobina_calc_y_cobertura():
    coil = E["propiedades"]["copper_coil_homog"]
    assert coil["estado"] == "calc_declarado"
    assert 0 < coil["k_W_mK"] < 390  # mezcla < Cu puro
    assert E["cobertura_inventario_f75"]["masa_cubierta_pct"] > 90
    assert any("[EST]" in a or "ANISOTROPO" in a for a in E["advertencias"])
    reg = EV / "registro_fem.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
