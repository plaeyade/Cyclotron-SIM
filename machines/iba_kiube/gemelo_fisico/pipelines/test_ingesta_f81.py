#!/usr/bin/env python3
"""Tests F81 — ingesta kiube v0.35 y verificación cruzada (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, os, pathlib

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072")) / "fase81"
E = json.load(open(EV / "ingesta_kiube_v035_f81.json"))
ETQ = "fuente_modelo_oficial_IBA (vía usuario)"


def test_integridad_md5_multiformato():
    assert E["fase"] == 81 and E["pq"] == "no_declarada_nunca"
    assert E["integridad"]["ok"] and E["integridad"]["verificados"] == E["paquete"]["n_artefactos"]
    f = E["paquete"]["por_formato"]
    for fmt in ("STEP AP214", "DAGMC h5m (MOAB)", "GDML (Geant4)", "STL binario"):
        assert f.get(fmt, 0) >= 1, fmt
    assert E["paquete"]["por_categoria"].get("telemetria_simulada", 0) >= 2
    assert E["paquete"]["por_categoria"].get("curva", 0) >= 5
    assert len(E["paquete"]["md5_manifiesto"]) == 32


def test_reetiquetado_honesto():
    v = E["verificacion"]
    # V04: radio de cámara, n.º, radio y posición de bocas coinciden con cotas de fuente IBA (T04F2/TUD03)
    assert set(v["re_etiquetadas"]) == {"V04-1", "V04-3", "V04-4", "V04-5"}
    reg = [json.loads(l) for l in (EV / "registro_reetiquetado.jsonl").read_text().splitlines() if l.strip()]
    for r in reg:
        if r["accion"] == "re_etiquetado":
            assert r["etiqueta_nueva"] == ETQ and r["etiqueta_anterior"].startswith("supuesto_declarado")
            assert r["valor_anterior"] is not None          # el supuesto anterior queda en el registro
        else:
            assert r["etiqueta_nueva"] == r["etiqueta_anterior"]   # nada se sobrescribe sin coincidencia + fuente IBA


def test_discrepancias_cuantificadas():
    d = {x["id"]: x for x in E["verificacion"]["discrepancias"]}
    assert -0.25 < d["V04-2"]["delta_rel"] < -0.15            # gap efectivo 132 → ~103 mm
    assert d["V04-7"]["delta_rel"] < -0.4                       # S_eff
    assert d["V04-11"]["delta_rel"] > 100                       # pérdida por gas: F65 frente a la medida de IBA
    assert 0.0 < d["F75-1"]["delta_rel"] < 0.05                 # acero 14.6 → 15.0 t
    assert "OM-3" in d                                          # amplificador 16 kW frente a la potencia del modelo
    for x in d.values():
        assert x["delta"] is not None and x["delta_rel"] is not None
    assert "F75-2" in E["verificacion"]["coinciden"]            # bobinas 2.12 t
    assert any(x["id"] == "V04-6" for x in E["verificacion"]["no_comparables"])


def test_impacto_en_cadena():
    i = E["impacto"]
    assert i["camara_difiere"] and i["voxel_afectado"]
    alta = {x["fase"] for x in i["afectadas"] if x["prioridad"] == "alta"}
    assert {"F66", "F74", "F70"} <= alta
    marcas = [json.loads(l) for l in (EV / "marcas_afectacion.jsonl").read_text().splitlines() if l.strip()]
    assert {m["fase_afectada"] for m in marcas} >= {"F66", "F74", "F70"}
    assert all(m["estado"] == "afectada_pendiente_re_corrida" for m in marcas)


def test_telemetria_simulada_declarada():
    t = E["telemetria_simulada"]
    assert t["origen"] == "simulado_declarado" and t["origen_por_fila_ok"]
    assert "fit_factor" in t["calibration_layer"]["uso_prohibido"]
    assert any("memoria 3" in a for a in E["advertencias"])
