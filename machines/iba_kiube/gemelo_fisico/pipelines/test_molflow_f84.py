#!/usr/bin/env python3
"""Tests F84 — Molflow+ oficial con la cámara v0.35 y correcciones de los motores de vacío."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, os, pathlib

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072")) / "fase84"
E = json.load(open(EV / "molflow_v035_f84.json"))


def test_validacion_analitica_ambos_motores():
    a = E["validacion_analitica"]
    assert abs(a["molflow_sobre_QS"] - 1) < 0.02 and abs(a["tpmc_sobre_QS"] - 1) < 0.02
    assert a["antes_de_corregir"]["tpmc_sobre_QS"] > 1.9 and a["antes_de_corregir"]["molflow_sobre_QS"] < 0.1


def test_molflow_vs_tpmc_corregido():
    for g in ("H2O", "H2"):
        v = E["gases"][g]
        assert 0.85 < v["ratio_medio"] < 1.15, g                       # mismo mapa con los dos motores
        assert v["fraccion_absorbida"] > 0.95, g                        # Molflow absorbe todo lo inyectado (estanco)
    assert E["motor"]["md5_binario"] and len(E["motor"]["md5_binario"]) == 32


def test_stripping_y_correcciones():
    s = E["stripping"]
    assert 0.8 < s["cociente_molflow_F82"] < 1.25
    corr = [json.loads(l) for l in (EV / "correcciones_historicas.jsonl").read_text().splitlines() if l.strip()]
    fases = {c["fase"] for c in corr}
    assert {"F66/F82", "F67", "F74"} <= fases
    res = [json.loads(l) for l in (EV / "resolucion_marcas_f81.jsonl").read_text().splitlines() if l.strip()]
    assert any(r["fase_afectada"] == "F67" and r["estado"] == "re_corrida_en_F84" for r in res)
    assert any("memoria 3" in a for a in E["advertencias"])
