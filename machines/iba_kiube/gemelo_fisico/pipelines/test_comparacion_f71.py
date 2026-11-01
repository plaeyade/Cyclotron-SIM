#!/usr/bin/env python3
"""Tests F71 — comparación triple de blindaje (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase71")
E = json.load(open(EV / "comparacion_triple_f71.json"))


def test_campos_canonicos():
    assert E["fase"] == 71
    assert E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert "RPOS-PLANO (BIND)" in E["fuentes"]


def test_tres_aproximaciones_presentes():
    c = E["comparacion"]
    assert set(c) == {"F52_csg_bind_sab", "F62_casco_dagmc", "F70_kiube_voxel_bunker"}
    for k, v in c.items():
        assert "fuente" in v and "geometria" in v and "metrica" in v


def test_valores_en_rango_fisico():
    c = E["comparacion"]
    assert 0 < c["F52_csg_bind_sab"]["valor"] < 1
    assert 0 < c["F62_casco_dagmc"]["fraccion_fuga"] < 1
    assert 0 < c["F70_kiube_voxel_bunker"]["fuera_bloque_rel"] < 0.05
    assert c["F70_kiube_voxel_bunker"]["full_estado"] == "completada"


def test_metricas_no_homogeneas_declaradas():
    assert any("NO homogéneas" in a for a in E["advertencias"])
    al = E["alineacion_muros_oeste_este"]
    assert al["muro_principal_espesor_m"] == 2.11
    assert "oeste" in al["criterio"] and "este" in al["criterio"]
    reg = EV / "registro_comparacion.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
