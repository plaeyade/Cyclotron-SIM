#!/usr/bin/env python3
"""Tests F85 — cierre de la cadena de blindaje F70–F78 (kiube v0.35, Geant4 v2)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, os, pathlib

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072")) / "fase85"
E = json.load(open(EV / "cierre_blindaje_f85.json"))


def test_canon():
    assert E["fase"] == 85 and E["tipo"] == "simulacion_teorica" and E["pq"] == "no_declarada_nunca"
    assert E["entradas"]["corridas"] >= 5 and len(E["entradas"]["md5_resultado"]) == 32
    assert any("memoria 3" in a for a in E["advertencias"])


def test_f70_f71_f78_exterior_directo():
    assert E["F70"]["cumple_criterio_F70"]
    assert 18 < E["F71"]["tvl_medio_cm"] < 40
    for k, m in E["F78"]["muros"].items():
        assert m["exterior_cara_uSv_h"] > 0 and m["exterior_cara_se"] < m["exterior_cara_uSv_h"], k   # error < 100 %
        assert 1e-12 < m["T"] < 1e-6, k
    assert E["F78"]["estado"].startswith("sustituida")


def test_f77_mallas_y_f72():
    m = E["F77"]
    assert m["finitas"] and m["no_negativas"] and m["voxeles"] > 500_000
    assert m["voxeles_con_dosis"] > 100 and m["voxeles_std_menor_50pct"] > 0
    assert set(E["F72"]) == {"Co60", "Eu152", "Cu64", "Ag110m"}
    assert all(v["A_sat_Bq"] > 0 for v in E["F72"].values())


def test_marcas_cerradas():
    res = [json.loads(l) for l in (EV / "resolucion_marcas_f81.jsonl").read_text().splitlines() if l.strip()]
    est = {r["fase_afectada"] for r in res}
    assert {"F70", "F71", "F72", "F76", "F77", "F78"} <= est
