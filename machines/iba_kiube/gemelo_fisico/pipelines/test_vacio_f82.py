#!/usr/bin/env python3
"""Tests F82 — re-corrida F66/F74/F65 con la cámara y el bombeo de kiube v0.35."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, os, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]
from app.domain.gas_stripping import MBAR_A_PA, densidad_gas_cm3            # noqa: E402

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072")) / "fase82"
E = json.load(open(EV / "vacio_v035_f82.json"))


def test_unidades_mbar_pa():
    assert MBAR_A_PA == 100.0
    # 1e-6 mbar a 293 K ≈ 2,5e10 moléculas/cm³ (valor de libro)
    assert 2.3e10 < densidad_gas_cm3(1e-6, 293.15) < 2.6e10


def test_entradas_v035_y_trazabilidad():
    e = E["entradas"]
    assert abs(e["gap_mm"] - 102.6) < 0.5 and e["r_cam"] == 548.0 and e["r_boca"] == 127.0
    assert all(len(h) == 32 for h in e["md5"].values())
    lados = sorted(b["lado"] for b in E["F66_tpmc"]["H2O"]["bombas"])
    assert lados == [-1, -1, 1, 1]                                   # 2 abajo (difusoras) + 2 arriba (carcasas RF)


def test_bombeo_analitico_coherente_con_kiube():
    for g in ("S_H2O_l_s", "S_H2_l_s"):
        f82, kiube = E["consistencia_con_kiube"][g]
        assert abs(f82 / kiube - 1) < 0.05


def test_tpmc_vs_analitico_y_medida():
    for g in ("H2O", "H2"):
        a = E["F74_analitico"][g]
        # F84: con el área de las dos tapas corregida el TPMC queda en ~1,0–1,1 × Q/S (antes ~2,0–2,2 por el error)
        assert 0.9 < a["tpmc_media_sobre_analitica"] < 3.0
        assert abs(E["F66_tpmc"][g]["captura_media_tpmc"] / a["s_captura"] - 1) < 0.05
    assert 0.4 < E["comparacion_V04"]["P_H2O_tpmc_sobre_medida_TUD03T"] < 1.5      # F84: 0,54 con el motor corregido


def test_stripping_frente_a_medida_iba():
    s = E["F65_stripping"]
    assert 0.02 < s["KIUBE180_180uA"]["fraccion_perdida"] < 0.25      # IBA: 19 % de pérdida total (TUD03); F84: 4,4 %
    assert s["sensibilidad_H2_3p5sccm_AMIT_pct"] > s["KIUBE180_180uA"]["pct"]
    c = E["comparacion_V04"]
    assert c["factor_por_unidades"] == 1000.0 and c["factor_por_geometria_bombeo_y_carga_v035"] > 1


def test_registros_append_only():
    corr = [json.loads(l) for l in (EV / "correcciones_historicas.jsonl").read_text().splitlines() if l.strip()]
    assert any(c["fase"] == "F65" and c["factor"] == 1000.0 for c in corr)
    res = [json.loads(l) for l in (EV / "resolucion_marcas_f81.jsonl").read_text().splitlines() if l.strip()]
    est = {r["fase_afectada"]: r["estado"] for r in res}
    assert est["F66"] == est["F74"] == est["F65"] == "re_corrida_en_F82"
    assert est["F67"] == "afectada_pendiente_re_corrida"
