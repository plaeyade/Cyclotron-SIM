#!/usr/bin/env python3
"""Tests F70 — sin numpy (stdlib), ejecutables con python3.11+pytest o python3.12."""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, struct

F70 = pathlib.Path(f"{_RPOS_RAIZ}/work/f70_mc")
EVI = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase70")

def leer_voxel(path):
    d = open(path, "rb").read()
    nx, ny = struct.unpack_from("<ii", d, 0)
    step, H, Z0, Z1 = struct.unpack_from("<dddd", d, 8)
    nm, = struct.unpack_from("<i", d, 40)
    nz = int((Z1 - Z0) / step)
    pos, counts = 44, {}
    while True:
        i, j = struct.unpack_from("<ii", d, pos)
        if not (0 <= i < nx and 0 <= j < ny):
            break
        buf = struct.unpack_from("<%dh" % nz, d, pos + 8)
        for v in buf:
            if v >= 0:
                counts[v] = counts.get(v, 0) + 1
        pos += 8 + 2 * nz
    names = []
    for _ in range(nm):
        L, = struct.unpack_from("<i", d, pos); pos += 4
        names.append(d[pos:pos + L].decode()); pos += L
    return {names[k]: v * step ** 3 for k, v in counts.items()}

def test_voxel_2cm_fracciones_volumen():
    vox = leer_voxel(F70 / "voxel_map_2cm.bin")
    vj = json.load(open(f"{_RPOS_RAIZ}/work/kiube_v08/kiube/out/mc/volumes.json"))
    ref = {}
    for info in vj.values():
        if isinstance(info, dict) and info.get("material") and info.get("volume_cm3"):
            ref[info["material"]] = ref.get(info["material"], 0.0) + info["volume_cm3"]
    assert 0.97 < vox["steel_1010"] / ref["steel_1010"] < 1.03      # yugo: estructural
    assert 0.95 < vox["stainless_304L"] / ref["stainless_304L"] < 1.05
    assert 0.93 < vox["copper_coil_homog"] / ref["copper_coil_homog"] < 1.10

def test_evidencia_canonica():
    e = json.load(open(EVI / "mc_maquina_bunker_f70.json"))
    assert e["fase"] == 70 and e["version"] == "2.6.0"
    assert e["tipo"] == "simulacion_teorica"
    assert e["pq"] == "no_declarada_nunca"
    for fuente in ("BON21", "RPOS-PLANO", "ICRP74", "ENDFB8", "KIUBE08"):
        assert fuente in e["fuentes"]

def test_verif_existe_y_es_honesta():
    v = json.load(open(F70 / "verif_voxel.json"))
    assert v["source_neutrons_per_s"] > 0
    assert 0.0 < v["dose_integral_fuera_bloque_rel"] < 0.05   # bunker contiene
    e = json.load(open(EVI / "mc_maquina_bunker_f70.json"))
    if e["corridas"]["estado_full"] == "pendiente_entorno":
        assert "pendiente" in e and e["pendiente"], "full pendiente debe declararse"

def test_registro_append_only():
    reg = EVI / "registro_mc.jsonl"
    assert reg.exists()
    lineas = [l for l in open(reg) if l.strip()]
    assert lineas, "registro vacio"
    r = json.loads(lineas[-1])
    assert r["fase"] == 70 and r["tipo"] == "simulacion_teorica"
    assert len(r["md5_voxel_2cm"]) == 32
