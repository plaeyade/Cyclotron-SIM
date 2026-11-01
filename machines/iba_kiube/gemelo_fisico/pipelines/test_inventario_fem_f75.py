#!/usr/bin/env python3
"""Tests F75 — inventario FEM (stdlib only)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase75")
E = json.load(open(EV / "inventario_fem_f75.json"))


def test_campos_canonicos():
    assert E["fase"] == 75 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert "INVENTARIO" in E["alcance"]


def test_inventario_completo_y_masas():
    vols = json.load(open(f"{_RPOS_RAIZ}/work/kiube_v08/kiube/out/mc/volumes.json"))
    assert E["n_componentes"] == len(vols) == len(E["componentes"])
    # masa de cada componente = V * rho (cuando hay densidad)
    mats = json.load(open(f"{_RPOS_RAIZ}/work/f70_mc/materials.json"))
    for nombre, c in E["componentes"].items():
        d = mats.get(c["material"], {}).get("density")
        if d:
            assert abs(c["masa_kg"] - c["volumen_cm3"] * d / 1000) < 1e-3
    # masa total = suma por material
    assert abs(E["masa_total_kg"] - sum(e["masa_kg"] for e in E["por_material"].values())) < 0.2
    # el concreto del búnker domina la masa total; el acero del yugo domina en la máquina
    top = max(E["por_material"].items(), key=lambda kv: kv[1]["masa_kg"])
    assert top[0] == "concrete_portland"
    assert E["por_material"]["steel_1010"]["masa_kg"] > 10000  # yugo > 10 t


def test_honestidad_propiedades_y_laminas():
    for c in E["componentes"].values():
        assert c["propiedades_fem"]["k_termica"] == "pendiente_fuente"
    assert E["clases_fem"].get("lamina_no_resuelta_declarada", 0) >= 1
    assert any("diferido a versión final" in p for p in E["pendientes"])
    assert any("BLOQUEADAS hasta haz" in p for p in E["pendientes"])


def test_caso_elmer_y_registro():
    assert E["caso_fem_existente"]["tipo"].startswith("Heat Equation")
    assert len(E["caso_fem_existente"]["md5"]) == 32
    assert E["cargas_termicas_con_fuente"]["perdida_haz_gas_w"]["valor"] == 0.0121
    reg = EV / "registro_fem.jsonl"
    assert reg.exists() and len(reg.read_text().strip().splitlines()) >= 1
