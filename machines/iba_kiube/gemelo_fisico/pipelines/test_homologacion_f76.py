#!/usr/bin/env python3
"""Tests F76 — homologación H*(10) (stdlib + numpy)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase76")
E = json.load(open(EV / "homologacion_h10_f76.json"))


def test_campos_canonicos():
    assert E["fase"] == 76 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert "ICRP-116" in E["metrica_homologada"]


def test_tres_geometrias_misma_metrica():
    t = E["tabla"]
    assert set(t) == {"F52_csg_bind_sab", "F62_casco_dagmc", "F70_kiube_voxel"}
    # F52: T_H10 coherente con T_flujo histórica (±30%, endurecimiento moderado)
    r = t["F52_csg_bind_sab"]["T_H10_muro_este"] / t["F52_csg_bind_sab"]["T_flujo_referencia_F52"]
    assert 0.7 < r < 1.3
    # F62: reproducción del flujo de muro original <0.5%
    f62 = t["F62_casco_dagmc"]
    assert abs(f62["flux_muro_media"] - 595.9) / 595.9 < 0.005
    assert f62["H10_muro_sigma_semillas"] / f62["H10_muro_media_pSv_cm3_part"] < 0.05


def test_h_efectivo_mismo_orden():
    h52 = E["tabla"]["F52_csg_bind_sab"]["h_efectivo_losa_interior_pSv_cm2"]
    h62 = E["tabla"]["F62_casco_dagmc"]["h_efectivo_muro_pSv_cm2"]
    assert h52 > 0 and h62 > 0
    assert 0.1 < h52 / h62 < 10  # mismo orden de magnitud


def test_cero_f70_declarado_y_registro():
    f70 = E["tabla"]["F70_kiube_voxel"]
    assert f70["H10_losa_exterior"] == 0.0
    assert "estadística insuficiente DECLARADA" in f70["nota_exterior_cero"]
    assert any("no cero físico" in a for a in E["advertencias"])
    reg = EV / "registro_homologacion.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
