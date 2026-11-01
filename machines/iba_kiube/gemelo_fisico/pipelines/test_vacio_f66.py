#!/usr/bin/env python3
"""Tests F66 — motor TPMC de vacío y acoplamiento con stripping."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia externa
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

from app.domain.vacuum_tpmc import CamaraVacio, simular_mapa_presion
from app.domain.gas_stripping import PerfilPresion, simular_perdida_gas

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase66")


def test_conservacion_todas_absorbidas():
    # en régimen molecular con bombas, toda molécula termina absorbida
    cam = CamaraVacio()
    r = simular_mapa_presion(cam, n_moleculas=20_000, semilla=1)
    assert r["absorbidos"] == 20_000
    # fracción de captura en boca ≈ s declarado
    assert abs(r["absorbidos"] / r["hits_bomba"] - cam.s_captura) < 0.02


def test_linealidad_en_q():
    # TPMC es lineal: duplicar Q duplica P (misma geometría/semilla)
    r1 = simular_mapa_presion(CamaraVacio(q_h2_fuentes_total_mbar_l_s=1e-3),
                              n_moleculas=20_000, semilla=7)
    r2 = simular_mapa_presion(CamaraVacio(q_h2_fuentes_total_mbar_l_s=2e-3),
                              n_moleculas=20_000, semilla=7)
    p1 = sum(r1["p_mbar"]) / len(r1["p_mbar"])
    p2 = sum(r2["p_mbar"]) / len(r2["p_mbar"])
    assert abs(p2 / p1 - 2.0) < 0.15


def test_gradiente_centro_a_pared():
    # fuente central + bombas periféricas → P máxima en el centro
    r = simular_mapa_presion(CamaraVacio(), n_moleculas=60_000, semilla=3)
    assert r["p_mbar"][0] > r["p_mbar"][-1]
    # mínimo cerca de la posición de las bombas (r=375 → bin ~13-14)
    imin = min(range(len(r["p_mbar"])), key=lambda i: r["p_mbar"][i])
    assert 11 <= imin <= 17


def test_p_mayor_que_q_s():
    # con conductancias finitas, P media > Q/S nominal
    r = simular_mapa_presion(CamaraVacio(), n_moleculas=30_000, semilla=5)
    p_media = sum(r["p_mbar"]) / len(r["p_mbar"])
    assert p_media > r["p_equilibrio_q_s_mbar"]


def test_perfil_desde_mapa_interpola():
    pf = PerfilPresion.desde_mapa([0.0, 100.0, 200.0], [1e-6, 5e-7, 2e-7])
    assert pf.presion_mbar(0.0) == 1e-6
    assert abs(pf.presion_mbar(150.0) - 3.5e-7) < 1e-12
    assert pf.presion_mbar(999.0) == 2e-7  # extrapolación plana


def test_perdida_con_mapa_coherente():
    pf = PerfilPresion.desde_mapa([0.0, 548.0], [4.5e-7, 1.4e-7])
    r = simular_perdida_gas(corriente_ua=200, energia_final_mev=18.0,
                            delta_e_kev_vuelta=98.5, b_media_t=1.360,
                            perfil=pf)
    # F82: corregida la conversión mbar → Pa (×1000); antes 1e-6 – 1e-3
    assert 1e-3 < r.fraccion_perdida < 1e-1


def test_f66_json_coherente():
    d = json.loads((EVID / "mapa_vacio_f66.json").read_text())
    assert "NO_produccion_real" in d["estado"]
    assert len(d["mapa_presion"]["p_mbar"]) == 20
    assert d["mapa_presion"]["n_moleculas"] == 600_000
    assert "TPMC-IPAC14" in d["fuentes"]["algoritmo"]
    assert any("Molflow" in p for p in d["pendientes"])
    res = d["resultados_perdida"]
    assert res["mapa_tpmc_f66"]["fraccion_perdida"] > 0


def test_registro_append_only():
    reg = EVID / "registro_vacio.jsonl"
    lineas = [l for l in reg.read_text().strip().splitlines() if l]
    assert len(lineas) >= 1
    e = json.loads(lineas[-1])
    assert e["fase"] == "F66" and e["p_media_mbar"] > 0
