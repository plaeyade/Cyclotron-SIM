#!/usr/bin/env python3
"""Tests F72 — activación y dosis al gemelo (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, math

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase72")
E = json.load(open(EV / "activacion_gemelo_f72.json"))


def test_campos_canonicos():
    assert E["fase"] == 72 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert "Informe RadiopharmaOS v2.0 Digital Twins" in E["fuentes"]


def test_isotopos_y_canales():
    iso = E["isotopos"]
    assert set(iso) == {"Co60", "Eu152", "Ag110m", "Cu64"}
    assert len(iso["Co60"]["canales"]) == 2  # acero (n,gamma) + cobre (n,a)
    for v in iso.values():
        assert v["A_sat_Bq"] > 0 and v["A_sat_Bq_sigma"] > 0
        assert set(v["A_EOB_Bq"]) == {"1h", "8h", "30d", "1a"}


def test_fisica_crecimiento():
    # A(t) monótona creciente y acotada por saturación
    for v in E["isotopos"].values():
        a = v["A_EOB_Bq"]
        assert a["1h"] <= a["8h"] <= a["30d"] <= a["1a"] <= v["A_sat_Bq"] * 1.0000001
        # vida corta: Cu64 a 30d debe estar casi en saturación
    cu = E["isotopos"]["Cu64"]
    assert cu["A_EOB_Bq"]["30d"] > 0.999 * cu["A_sat_Bq"]


def test_honestidad_y_registro():
    assert E["dosis_agregada"]["estado_mallas_por_voxel"].startswith("pendiente_entorno")
    assert any("TEÓRICO" in a or "teórico" in a.lower() for a in E["advertencias"])
    assert any("UNIDADES PENDIENTE" in a for a in E["advertencias"])
    assert E["verificacion_unidades"].startswith("pendiente")
    assert abs(E["fuente_neutrones_n_s"] - 4616232340295.063) < 1
    reg = EV / "registro_activacion.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
