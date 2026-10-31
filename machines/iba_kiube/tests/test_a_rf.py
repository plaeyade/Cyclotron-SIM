"""Pruebas de RF (lentas). Van en un archivo propio que pytest corre PRIMERO, antes de que las pruebas de geometría
carguen cientos de piezas CAD en memoria (con 7 GB, juntas no caben)."""
import json, subprocess, sys
import pytest
from conftest import ROOT


@pytest.mark.lento
def test_rf_soluciones_exactas():
    r = subprocess.run([sys.executable, str(ROOT / "fem/rf_eigen.py"), "--verificar"], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-1500:]
    for k, v in json.load(open(ROOT / "validation/rf_verificacion.json")).items():
        assert abs(v["error_rel"]) < 1e-4, (k, v)


@pytest.mark.lento
def test_rf_frecuencia_del_modo_acelerador(PE):
    P, _, _ = PE
    r = subprocess.run([sys.executable, str(ROOT / "fem/rf_kiube.py"), "--maxh", "0.035", "--tag", "prueba"],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-1500:]
    R = json.load(open(ROOT / "out/fem_rf/resultado_prueba.json"))
    m0 = R["modos"][0]; f = m0["f_MHz"]
    assert m0["fase"] == "en fase" and abs(m0["V1_sobre_V0"] - 1) < 0.05, m0       # dees equilibrados
    assert abs(f / P["beam"]["rf_MHz"] - 1) < 0.005, f"f = {f:.3f} MHz"
    if P["dee"].get("puente_w"):                                                   # v0.19: dees acoplados en el centro
        otros = [m["f_MHz"] for m in R["modos"][1:]]
        assert all(abs(x - f) > 10 for x in otros), otros                          # modo acelerador aislado
    if "acoplador" in R:                                                            # lazo crítico de tamaño razonable
        assert 5 < R["acoplador"]["lazo_h_critico_mm"] < 150, R["acoplador"]
