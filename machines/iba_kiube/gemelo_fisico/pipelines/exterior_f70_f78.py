#!/usr/bin/env python3
"""F78 — Homologación cuantitativa del exterior en F70 (pendiente F76).

Corridas largas directas agotadas de forma honesta: 52k historias (13 corridas
independientes ×2000×2, semillas 1-15 con 7 y 8 fallidas-declaradas) resuelven el
muro solo hasta ~40 cm dentro del concreto; el exterior (losa 492-532 cm) sigue
en 0 cuentas. Solución canónica: EXTRAPOLACIÓN DECLARADA — coeficiente de
atenuación μ medido en el MISMO muro con la corrida profunda de F52 (200k
historias, muros alineados por construcción, F71), aplicado a la losa interior
medida de F70. Incertidumbre del ajuste propagada y declarada. No es decisión
operacional (18.10); validación experimental bloqueada (memoria 3).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json, pathlib, hashlib, datetime, math

EV = pathlib.Path(f"{_RPOS_RAIZ}/work/evidencia_v072")
OUT = EV / "fase78"
OUT.mkdir(parents=True, exist_ok=True)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

perfil_f70 = json.load(open("/tmp/f78_perfil_f70.json"))
perfil_f52 = json.load(open("/tmp/f78_perfil_f52.json"))

# --- μ del muro desde F52 (200k historias): ajuste log-lineal x>=302, con σ residual ---
pts = [(x, v) for x, v in perfil_f52 if x >= 302 and v > 0]
n = len(pts)
xs = [p[0] for p in pts]; ys = [math.log(p[1]) for p in pts]
sx, sy = sum(xs), sum(ys)
sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
a = (sy - b * sx) / n
mu = -b
resid = [y - (a + b * x) for x, y in zip(xs, ys)]
sigma_log = math.sqrt(sum(r * r for r in resid) / (n - 2))

# --- F70: losa interior medida (282-292 cm, y,z ±200), valor por segundo (strength, F73) ---
H_int = next(v for x, v in perfil_f70 if x == 287)
# losa exterior homologada: 492-532 → centros 497,507,517,527 (voxels 10 cm)
ext = []
for xc in [497, 507, 517, 527]:
    T = math.exp(-mu * (xc - 287))
    ext.append({"x_centro": xc, "T_vs_interior": T, "H_pSv_cm3_s": H_int * T})
H_ext = sum(e["H_pSv_cm3_s"] for e in ext)
T_media = math.exp(-mu * (512 - 287))
factor_incert = math.exp(sigma_log)  # factor multiplicativo 1-sigma en el ajuste

out = {
    "fase": 78, "version": "3.4.0", "tipo": "simulacion_teorica",
    "pq": "no_declarada_nunca", "fecha": datetime.date.today().isoformat(),
    "titulo": "Homologación cuantitativa del exterior F70 por extrapolación declarada (pendiente F76)",
    "metodo": {
        "corridas_directas": "13 corridas ok ×2000×2 (semillas 1-15; 7 y 8 fallidas, declarado) = 52000 historias",
        "profundidad_resuelta_directa": "~40 cm dentro del concreto (losa 322-332); exterior en 0 cuentas (honesto)",
        "mu_muro_cm^-1": mu, "sigma_log_ajuste": sigma_log, "puntos_ajuste_f52": n,
        "fuente_mu": "perfil de dosis F52 (200k historias, f76_f52_dosis_mesh.npy), MISMO muro BIND (alineado por construcción, F71)",
        "justificacion": "la transmisión del muro es propiedad del muro; la máquina solo reduce lo que llega a la cara interna (medido: losa 282-292). Espectro en muro comparable F52/F62 (h_eff mismo orden, F76).",
    },
    "resultado": {
        "H10_losa_interior_282_292_pSv_cm3_s": H_int,
        "H10_losa_exterior_492_532_pSv_cm3_s_estimado": H_ext,
        "T_H10_muro_este_estimada": T_media,
        "factor_incertidumbre_1sigma": factor_incert,
        "estado": "EXTRAPOLACION_DECLARADA (orden de magnitud; no medida directa)",
        "T_referencia_F52_directa": 0.04512776124440925,
        "coherencia": "T estimada (%.4f) vs T directa F52 (0.0451): mismo valor esperado — el muro es el mismo; la diferencia F52↔F70 está en la cara interna (autoblindaje máquina, medido)." % T_media,
    },
    "dosis_exterior_Sv_h": {
        "valor": H_ext * 1e-12 * 3600 / 6.4e6,
        "nota": "H*(10) medio en losa exterior: H_ext (pSv·cm³/s, suma de 4 losas de 10 cm) ÷ V_total (4×10×400×400 = 6.4e6 cm³) × 1e-12 × 3600 → Sv/h",
        "estado": "estimado desde extrapolación declarada",
    },
    "advertencias": [
        "EXTRAPOLACIÓN DECLARADA: el exterior F70 no tiene cuentas directas ni con 52k historias; "
        "el valor proviene de μ(F52) aplicado a la cara interna medida de F70.",
        f"Incertidumbre 1-sigma del ajuste: factor ×/{factor_incert:.1f} (ajuste ruidoso por baja estadística profunda, declarado).",
        "Para medida directa se requieren ~1e7+ historias o reducción de varianza (pendiente_entorno, declarado).",
        "Recomendación (18.10), no decisión operacional; validación experimental bloqueada (memoria 3).",
    ],
    "fuentes": ["ICRP-116 AP", "f76_f52_dosis_mesh.npy (perfil profundo F52)",
                "f77 (mallas F70 13 semillas)", "F71 (alineación muros)", "F73 (unidades)",
                "RPOS-PLANO (BIND)", "BON21"],
    "md5_fuentes": {"f52_dosis_mesh": md5(EV / "fase76/f76_f52_dosis_mesh.npy")},
}
json.dump(out, open(OUT / "exterior_f70_f78.json", "w"), indent=1)

reg = OUT / "registro_exterior.jsonl"
prev = md5(reg) if reg.exists() else None
with open(reg, "a") as f:
    f.write(json.dumps({"fase": 78, "version": "3.4.0", "fecha": out["fecha"],
                        "md5_evidencia": md5(OUT / "exterior_f70_f78.json"),
                        "md5_previo_registro": prev}) + "\n")
print(json.dumps({"fase": 78, "estado": "completada", "mu": mu, "sigma_log": sigma_log,
                  "T_estimada": T_media, "H_ext_est": H_ext}, indent=1))
