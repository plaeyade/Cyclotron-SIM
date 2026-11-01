#!/usr/bin/env python3
"""F70 — MC máquina KIUBE voxelizada (2 cm) + búnker BIND (CSG). OpenMC.

Cierra la fase con lo verificado y declara lo pendiente. Nunca fabrica
resultados: la corrida full queda pendiente_entorno (inestabilidad del
sandbox, 2026-10-18) con reanudación automática documentada.

Canon: simulacion_teorica; fuentes citadas (49.11); PQ no_declarada_nunca
(14.9); registro append-only ALCOA+.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, struct, subprocess, sys

HERE = pathlib.Path(__file__).resolve().parent          # repo/pipelines
F70 = pathlib.Path(f"{_RPOS_RAIZ}/work/f70_mc")
EVI = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase70")
EVI.mkdir(parents=True, exist_ok=True)

def md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()

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
    return {names[k]: v * step ** 3 for k, v in counts.items()}, (nx, ny, nz, step)

def ejecutar():
    # 1) validación de fracciones de volumen (voxel 2 cm vs CAD volumes.json)
    vj = json.load(open(f"{_RPOS_RAIZ}/work/kiube_v08/kiube/out/mc/volumes.json"))
    ref = {}
    for info in vj.values():
        if isinstance(info, dict) and info.get("material") and info.get("volume_cm3"):
            ref[info["material"]] = ref.get(info["material"], 0.0) + info["volume_cm3"]
    vox2, dim2 = leer_voxel(F70 / "voxel_map_2cm.bin")
    vox1, dim1 = leer_voxel(F70 / "voxel_map.bin")
    comp = {}
    for m in ("steel_1010", "copper_coil_homog", "stainless_304L", "copper_ofe", "niobium"):
        comp[m] = {"cad_cm3": ref[m], "vox2cm_cm3": vox2.get(m, 0.0),
                   "ratio_2cm": vox2.get(m, 0.0) / ref[m],
                   "ratio_1cm": vox1.get(m, 0.0) / ref[m]}
    no_resueltos = ["havar", "titanium", "graphite_foil", "water_o18(parcial)", "helium(parcial)"]

    # 2) resultado de verificación (única corrida completada y guardada)
    verif = json.load(open(F70 / "verif_voxel.json"))

    full = F70 / "full_voxel.json"
    estado_full = "completada" if full.exists() else "pendiente_entorno"
    full_res = json.load(open(full)) if full.exists() else None

    out = {
        "fase": 70, "version": "2.6.0",
        "tipo": "simulacion_teorica",
        "pq": "no_declarada_nunca",
        "titulo": "MC maquina KIUBE (voxel 2 cm) + bunker BIND (CSG), OpenMC",
        "geometria": {
            "ruta": "voxelizacion MOAB OBB ray-casting + CSG",
            "motivo_voxel": "61 solapes reales en kiube_openmc.h5m (overlap_check DAGMC); el '0 solapes' de kiube era muestreo G4 de 2000 pts/vol, insuficiente",
            "voxel_cm": 2.0,
            "celosia": {"dims": dim2, "margen_vacio": "1 voxel (bug segfault borde)"},
            "solapes_resueltos_por": "prioridad volumen menor (declarado)",
            "bunker_csg": "sala 5.9x4.7x2.5 m; bloque hormigon +/-506/+/-446/-167..283 cm; muro 2.11 m",
            "no_resueltos_2cm": no_resueltos,
            "validacion_volumenes": comp,
        },
        "fuentes": {
            "BON21": "Bonvin 2021 Fig.7 — espectro n H2-18O (238 ptos), 4.93e-3 n/p",
            "RPOS-PLANO": "plano BIND — geometria bunker",
            "ICRP74": "ICRP-74/116 — coeficientes H*(10) AP",
            "ENDFB8": "ENDF/B-VIII.0 — 67 nucleidos (solo los del modelo)",
            "KIUBE08": "kiube v0.8 — modelo, materiales, fuentes, tallies",
        },
        "corridas": {
            "verificacion_1cm_3000x2": verif,
            f"full_{full_res.get("corridas_ok","?")}x2000x2_malla10cm": full_res,
            "estado_full": estado_full,
            "nota_estabilidad": "openmc 0.16.0 segfault >~12000 historias/corrida en esta celosia; 10k x 2 estable (verificado experimentalmente)",
            "reanudacion": "sh /mnt/agents/work/f70_mc/f70_master.sh + f70_full.py (reanudable, progreso en estado/f70_progreso.npz)",
        },
        "advertencias": [
            "verificacion de baja estadistica (3000x2): valores de dosis/reaccion NO usar cuantitativamente hasta completar la full",
            "estadisticas de la full (cuando exista): 20 corridas independientes 10k x 2, media de medias",
        ],
        "pendiente": [] if full.exists() else [
            "corrida full (pendiente_entorno: reinicios del sandbox 2026-10-18; reanudar >= 2026-11-04)",
            "F71: comparacion triple kiube-voxel vs F62 casco vs F52 CSG (alinear muro oeste/este)",
        ],
    }
    json.dump(out, open(EVI / "mc_maquina_bunker_f70.json", "w"), indent=1)

    # registro append-only ALCOA+
    reg = EVI / "registro_mc.jsonl"
    with open(reg, "a") as f:
        f.write(json.dumps({"fase": 70, "version": "2.6.0",
                            "md5_voxel_1cm": md5(F70 / "voxel_map.bin"),
                            "md5_voxel_2cm": md5(F70 / "voxel_map_2cm.bin"),
                            "estado_full": estado_full,
                            "tipo": "simulacion_teorica"}) + "\n")
    print(json.dumps({"fase": 70, "estado_full": estado_full,
                      "verif_fuera_bloque_rel": verif.get("dose_integral_fuera_bloque_rel")}, indent=1))
    return out

if __name__ == "__main__":
    ejecutar()
