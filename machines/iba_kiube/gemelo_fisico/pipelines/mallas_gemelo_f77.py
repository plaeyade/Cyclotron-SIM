#!/usr/bin/env python3
"""F77 — Mallas de dosis H*(10) y flujo por voxel para el gemelo (pendiente F72).

Regenera y persiste las mallas por voxel (10 cm) de la geometría F70, con media
y std entre 5 corridas independientes, empaquetadas (npz + metadatos JSON) para
el Physics Twin / Digital Facility Twin. Resuelve el pendiente_entorno de F72.
Reglas canónicas: simulacion_teorica, pq no_declarada_nunca, unidades declaradas
(normalización strength, F73), fuentes citadas.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime
import numpy as np

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase77"

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

z = np.load(OUT / "f77_mallas_gemelo.npz")
meta = json.load(open(OUT / "f77_mallas_gemelo_meta.json"))

dose_m, dose_s = z["dose_h10_media"], z["dose_h10_std"]
flux_m, flux_s = z["flux3g_media"], z["flux3g_std"]

resumen = {
    "dose_h10": {
        "voxels_con_dosis": int((dose_m > 0).sum()), "n_voxels": int(dose_m.size),
        "max_pSv_cm3_s": float(dose_m.max()), "suma_total": float(dose_m.sum()),
        "voxels_std_menor_50pct": int(((dose_s < 0.5 * dose_m) & (dose_m > 0)).sum()),
    },
    "flux3g": {
        "voxels_con_flujo": int((flux_m > 0).sum()), "n_bins": int(flux_m.size),
        "max_part_cm2_s": float(flux_m.max()), "suma_total": float(flux_m.sum()),
    },
}

out = {**meta,
       "resumen": resumen,
       "artefactos": {"npz": "f77_mallas_gemelo.npz (4 arrays: dose_h10_media/std, flux3g_media/std)",
                      "meta": "f77_mallas_gemelo_meta.json"},
       "estado_pendiente_f72": "RESUELTO (mallas por voxel persistidas en /mnt)"}
json.dump(out, open(OUT / "mallas_gemelo_f77.json", "w"), indent=1)

reg = OUT / "registro_mallas.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 77, "version": "3.3.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "mallas_gemelo_f77.json"),
                        "md5_npz": md5(OUT / "f77_mallas_gemelo.npz"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 77, "estado": "completada",
                  "dose_max": resumen["dose_h10"]["max_pSv_cm3_s"],
                  "voxels_con_dosis": resumen["dose_h10"]["voxels_con_dosis"]}, indent=1))
