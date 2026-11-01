#!/usr/bin/env python3
"""Tests F77 — mallas por voxel para el gemelo (stdlib + numpy)."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib
import numpy as np

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase77")
E = json.load(open(EV / "mallas_gemelo_f77.json"))
Z = np.load(EV / "f77_mallas_gemelo.npz")


def test_campos_canonicos_y_pendiente():
    assert E["fase"] == 77 and E["tipo"] == "simulacion_teorica"
    assert E["pq"] == "no_declarada_nunca"
    assert E["estado_pendiente_f72"].startswith("RESUELTO")
    assert "memoria 3" in E["advertencias"][0]


def test_arrays_completos_y_forma():
    assert set(Z.files) == {"dose_h10_media", "dose_h10_std", "flux3g_media", "flux3g_std"}
    n = E["malla"]["n_voxels"]
    assert Z["dose_h10_media"].size == n == 660000
    assert Z["flux3g_media"].size == 3 * n
    for k in Z.files:
        assert np.all(np.isfinite(Z[k])) and np.all(Z[k] >= 0)


def test_estadistica_declarada():
    r = E["resumen"]
    assert r["dose_h10"]["voxels_con_dosis"] > 100
    assert r["dose_h10"]["max_pSv_cm3_s"] > 0
    # std entre semillas no negativa y finita; algunos voxels con buena estadística
    assert r["dose_h10"]["voxels_std_menor_50pct"] > 0


def test_coherencia_con_f76():
    # la media de dosis debe coincidir con la malla guardada en F76 (mismas semillas/protocolo)
    f76 = np.load(f"{_RPOS_RAIZ}/work/evidencia_v072/fase76/f76_f70_dosis_mesh_media.npy")
    d = Z["dose_h10_media"]
    assert d.shape == f76.shape
    s1, s2 = d.sum(), f76.sum()
    assert abs(s1 - s2) / s2 < 0.05  # mismas semillas: casi idéntica
    reg = EV / "registro_mallas.jsonl"
    assert reg.exists() and "md5_npz" in reg.read_text()
