#!/usr/bin/env python3
"""Tests F83 — re-corrida de F70 con la geometría kiube v0.35 (Geant4)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, os, pathlib

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072")) / "fase83"
E = json.load(open(EV / "blindaje_v035_f83.json"))


def test_canon_y_trazabilidad():
    assert E["fase"] == 83 and E["tipo"] == "simulacion_teorica" and E["pq"] == "no_declarada_nunca"
    assert len(E["entradas"]["md5_mapa"]) == 32 and len(E["entradas"]["md5_resultado"]) == 32
    assert E["entradas"]["corridas"] >= 5                           # media y error entre corridas independientes
    assert E["entradas"]["fuente"]["puertas_HDPE"] is True
    assert any("memoria 3" in a for a in E["advertencias"])


def test_muros_fisicos():
    m = E["resultados"]["muros"]
    assert set(m) == {"oeste", "este", "sur", "norte"}
    for v in m.values():
        assert v["exterior_uSv_h"] > 0                              # cuentas directas fuera del muro (F70: 0)
        assert 1e-12 < v["T_directa"] < 1e-6                        # 2,11 m de hormigón
    assert 18 < E["resultados"]["tvl_medio_cm"] < 40                # SRG Tabla 9: TVL ≈ 30 cm


def test_contraste_F70_F78():
    c = E["comparacion"]
    assert c["F70_fuera_bloque_rel"]["cumple"]
    assert c["F78_extrapolacion"]["cociente_T_F78_sobre_F83"] > 1e3     # F78 sobreestima la transmisión
    assert "discrepancia" in c["F78_extrapolacion"]["veredicto"]


def test_activacion_y_marcas():
    a = E["resultados"]["activacion"]
    assert set(a) == {"Co60", "Eu152", "Cu64", "Ag110m"}
    for iso in ("Co60", "Eu152", "Cu64"):
        assert 0.1 < a[iso]["cociente_F83_F72"] < 10                # mismo orden que F72
    res = [json.loads(l) for l in (EV / "resolucion_marcas_f81.jsonl").read_text().splitlines() if l.strip()]
    est = {r["fase_afectada"]: r["estado"] for r in res}
    assert est["F70"] == "re_corrida_en_F83" and est["F78"] == "contradicha_por_F83"
