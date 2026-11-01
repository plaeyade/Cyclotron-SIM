"""Pruebas Fase 67 — validación cruzada Molflow oficial vs TPMC propio.

Simulación teórica (memoria 3). PQ no_declarada_nunca (Informe 14.9).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "pipelines"))
EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase67")
JSON67 = EVID / "molflow_xval_f67.json"
CSV = "/tmp/mfrun/out3/facet_physics.csv"
XML = "/tmp/kiube_camara.xml"


def test_generador_xml_existe_y_valido():
    from tools.kiube_camara_molflow import generar_xml  # noqa
    out = generar_xml("/tmp/test_camara_f67.xml")
    assert out["facetas"] > 1000 and out["bombas"] == 4
    root = ET.parse("/tmp/test_camara_f67.xml").getroot()
    assert root.tag == "SimulationEnvironment"
    assert root.find("Geometry/Structures").get("nb") == "1"


def test_estructura_xml_camara():
    root = ET.parse(XML).getroot()
    facets = root.find("Geometry/Facets")
    assert int(facets.get("nb")) == len(list(facets))
    sticks = [float(f.find("Sticking").get("constValue")) for f in facets]
    assert any(s > 0.05 for s in sticks), "debe haber puertos con sticking"
    for f in facets:
        st = f.find("Structure")
        assert st.get("inStructure") == "0" and st.get("linksTo") == "0"


def test_molflow_csv_coherente():
    import csv
    with open(CSV) as fh:
        filas = list(csv.DictReader(fh))
    assert len(filas) > 1000
    for r in filas[:50]:
        assert float(r["Pressure [mbar]"]) > 0


@pytest.mark.skipif(not JSON67.exists(), reason="pipeline aún no ejecutado")
def test_json_f67_campos():
    d = json.loads(JSON67.read_text())
    assert d["fase"] == 67 and d["estado"] == "simulacion_teorica"
    assert d["pq"] == "no_declarada_nunca"
    assert len(d["mapa_molflow_f67"]["p_mbar"]) >= 15
    assert d["fuentes"], "toda figura requiere fuente (49.11)"


@pytest.mark.skipif(not JSON67.exists(), reason="pipeline aún no ejecutado")
def test_gradiente_y_orden_magnitud():
    d = json.loads(JSON67.read_text())
    p = d["mapa_molflow_f67"]["p_mbar"]
    assert p[0] > p[-1], "gradiente centro>pared esperado"
    assert all(1e-9 < x < 1e-5 for x in p), "orden de magnitud UHV"


@pytest.mark.skipif(not JSON67.exists(), reason="pipeline aún no ejecutado")
def test_xval_razonable_mismo_orden():
    d = json.loads(JSON67.read_text())
    ratio = d["comparacion"]["ratio_medio_molflow_vs_propio"]
    assert 0.1 < ratio < 10.0, "mismo orden de magnitud entre motores"


@pytest.mark.skipif(not JSON67.exists(), reason="pipeline aún no ejecutado")
def test_stripping_con_mapa_molflow_positivo_y_acotado():
    d = json.loads(JSON67.read_text())
    for nombre, c in d["comparacion"]["perdida_stripping"].items():
        assert 0 < c["fraccion_perdida"] < 0.01, nombre
        assert c["estado"] == "simulacion_teorica"


@pytest.mark.skipif(not JSON67.exists(), reason="pipeline aún no ejecutado")
def test_registro_append_only():
    reg = EVID / "registro_vacio.jsonl"
    lineas = reg.read_text().strip().splitlines()
    assert len(lineas) >= 1
    for ln in lineas:
        assert json.loads(ln)["estado"] == "simulacion_teorica"
