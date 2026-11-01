#!/usr/bin/env python3
"""F85 — Cierre de la cadena de blindaje F70–F78 con la geometría kiube v0.35.

Toma la corrida de kiube de F83 ampliada con corridas v1b (mismo esquema de importancia que F83 + malla de fluencia de
neutrones en 3 grupos; la dosis combina todas las corridas, las mallas de fluencia salen de las v1b) y cierra cada fase
afectada. Un esquema v2 (importancia también en techo y losa) se probó y se descartó: no completaba 500 historias en
25 min (registro en kiube mc_benchmark/bunker_mapa_v2_descartado.cc). El artefacto de varianza de F83 sigue declarado.
  F70  MC máquina + búnker                → sustituida por F83/F85 (Geant4, geometría v0.35)
  F71  comparación de métricas de muro    → tabla F83 / F85 / SRG (Tabla 9) / benchmark SRG de kiube v0.27
  F72  activación                         → actualizada con error estándar entre corridas
  F76  homologación de la métrica         → H*(10) ICRP-74 (n + γ) en todas las cifras de F83/F85; E(AP) ICRP-116 no
                                            reproducida (coeficientes pendiente_fuente) — declarado
  F77  mallas para el gemelo              → mallas_gemelo_f77.npz (dose_h10_media/std, flux3g_media/std), mismo formato
  F78  exterior                           → medida DIRECTA (sin extrapolación) y dosis anual con la carga de trabajo
Reglas: simulacion_teorica, pq no_declarada_nunca, append-only, 18.10, memoria 3.
Rutas: KIUBE_RAIZ, RPOS_EVIDENCIA, F85_CASO (por defecto out/mc_f70v035/h2500_t1450).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", f"{_RPOS_RAIZ}")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", (_RPOS_RAIZ + ""))  # raíz de evidencia/datos externos
import datetime, hashlib, json, math, os, pathlib
import numpy as np

KIUBE = pathlib.Path(os.environ.get("KIUBE_RAIZ", (_RPOS_RAIZ + "/work/kiube_v035/kiube")))
CASO = KIUBE / os.environ.get("F85_CASO", "out/mc_f70v035/h2500_t1450")
EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", (_RPOS_RAIZ + "/work/evidencia_v072")))
OUT = EV / "fase85"
# carga de trabajo anual del proyecto (kiube v0.28, estudio del techo): 50–150 k µAh/año; KIUBE 180 → 180 µA en blanco
CARGA_uAh = (50e3, 150e3)
LIMITES = {"publico_mSv_a": 1.0, "zona_radiacion_uSv_h": 50.0}   # RGSR México Art. 38 y Art. 113 (kiube v0.28)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def ejecutar():
    OUT.mkdir(parents=True, exist_ok=True)
    R2 = json.load(open(CASO / "resultado_mapa.json"))
    f83 = json.load(open(EV / "fase83/blindaje_v035_f83.json")) if (EV / "fase83/blindaje_v035_f83.json").exists() else None
    Z = np.load(CASO / "mallas_gemelo_f77.npz")
    dose, dstd, flu, fstd = (Z[k] for k in ("dose_h10_media", "dose_h10_std", "flux3g_media", "flux3g_std"))
    mallas = dict(voxeles=int(dose.size), dims=list(dose.shape), finitas=bool(np.isfinite(dose).all() and np.isfinite(flu).all()),
                  no_negativas=bool((dose >= 0).all() and (flu >= 0).all()), voxeles_con_dosis=int((dose > 0).sum()),
                  voxeles_std_menor_50pct=int(((dose > 0) & (dstd < 0.5 * dose)).sum()),
                  max_uSv_h=float(dose.max()), md5=md5(CASO / "mallas_gemelo_f77.npz"))
    horas = [c / 180.0 for c in CARGA_uAh]
    muros, tabla_f71 = {}, []
    for k, c in R2["caras_plano_medio"].items():
        ext, se = c["exterior_0_30cm_cara_uSv_h"]; pe, pse = c["exterior_0_30cm_parche_uSv_h"]
        anual = [max(ext, pe) * h / 1000 for h in horas]               # mSv/año con el máximo (cara o parche)
        muros[k] = dict(exterior_cara_uSv_h=ext, exterior_cara_se=se, exterior_parche_uSv_h=pe, exterior_parche_se=pse,
                        interior_uSv_h=c["interior_cara_uSv_h"], T=ext / c["interior_cara_uSv_h"] if c["interior_cara_uSv_h"] else None,
                        tvl_cm=c["tvl_cara_cm"], dosis_anual_mSv=anual,
                        cumple_publico=all(a < LIMITES["publico_mSv_a"] for a in anual),
                        fuera_de_zona=max(ext, pe) < LIMITES["zona_radiacion_uSv_h"])
        v1 = f83["resultados"]["muros"][k] if f83 else None
        tabla_f71.append(dict(muro=k, F83_12corridas_exterior=v1 and v1["exterior_uSv_h"], F83_se=v1 and v1["exterior_se"],
                              F85_exterior=ext, F85_se=se, F83_tvl=v1 and v1["tvl_cm"], F85_tvl=c["tvl_cara_cm"],
                              coherente_F83_F85=(abs(ext - v1["exterior_uSv_h"]) < 3 * math.hypot(se, v1["exterior_se"])) if v1 else None))
    tvls = [m["tvl_cm"] for m in muros.values() if m["tvl_cm"]]
    act = R2["activacion_saturacion"]
    out = {
        "fase": 85, "version": "3.11.0", "tipo": "simulacion_teorica", "pq": "no_declarada_nunca",
        "fecha": datetime.datetime.now().isoformat(timespec="seconds"),
        "titulo": "Cierre de la cadena de blindaje F70–F78 con kiube v0.35 (Geant4, F83 ampliada)",
        "entradas": {"caso": str(CASO.relative_to(KIUBE)), "md5_resultado": md5(CASO / "resultado_mapa.json"),
                     "corridas": R2["corridas"], "historias": R2["historias_totales"], "fuente": R2["fuente"],
                     "bloque": R2["bloque_hormigon"]},
        "F70": {"estado": "sustituida_por_F83_F85", "fuera_bloque_rel": R2["integral_dosis_fuera_bloque_rel"],
                "cumple_criterio_F70": R2["integral_dosis_fuera_bloque_rel"] < 0.05},
        "F71": {"tabla_muros": tabla_f71, "nota_coherencia": "F85 contiene las corridas de F83: la coherencia no es independiente", "tvl_medio_cm": float(np.mean(tvls)) if tvls else None,
                "referencias": {"SRG_tabla9_tvl_cm": "≈ 30 (200 → 220 cm: ×4,6)", "benchmark_kiube_v027_tvl_cm": 28.6},
                "nota": "F52 (CSG) y F62 (DAGMC) no se recalculan aquí: su evidencia no está en este entorno"},
        "F72": {iso: dict(A_sat_Bq=v["A_sat_Bq"], se_Bq=v["sigma_Bq"], fuente_xs=v["fuente_xs"]) for iso, v in act.items()},
        "F76": {"metrica": "H*(10), coeficientes ICRP-74 (n y γ), en todas las cifras de F83/F85",
                "pendiente": "E(AP) ICRP-116: coeficientes pendiente_fuente; no se reproduce la homologación de F76"},
        "F77": {"archivo": str((CASO / "mallas_gemelo_f77.npz").relative_to(KIUBE)), **mallas,
                "unidades": ["µSv/h (KIUBE 180, 2 × 90 µA)", "n/cm²/s en 3 grupos"]},
        "F78": {"estado": "sustituida: exterior DIRECTO", "muros": muros,
                "techo_exterior_uSv_h": R2["techo_exterior_medio_uSv_h"],
                "carga_trabajo_uAh_a": CARGA_uAh, "horas_haz_a": horas, "limites": LIMITES,
                "cumple_todos_los_muros": all(m["cumple_publico"] and m["fuera_de_zona"] for m in muros.values())},
        "advertencias": [
            "Simulación teórica: sin haz real (memoria 3). No cierra 18.19; PQ no_declarada_nunca (14.9).",
            "Artefacto de varianza de F83 (meseta a 1,4–2,0 m) sigue presente: μ se ajusta en 40–130 cm y el exterior se mide 0–30 cm fuera.",
            "Techo 1,45 m y losa 0,55 m de F70 (comparabilidad); el techo real sigue pendiente de los cortes (SRG mín. 1,7 m).",
            "Dosis anual = tasa exterior máxima (cara o parche) × horas de haz; sin factor de ocupación (conservador).",
            "Geometría del búnker estimada sobre el plano BIND (muro oeste medido; resto EST igual).",
        ],
    }
    json.dump(out, open(OUT / "cierre_blindaje_f85.json", "w"), indent=1, ensure_ascii=False, default=float)
    with open(OUT / "resolucion_marcas_f81.jsonl", "a") as f:
        for fase, est in (("F70", "sustituida_por_F83_F85"), ("F71", "cerrada_en_F85 (tabla F83/F85/SRG)"),
                          ("F72", "actualizada_en_F85"), ("F76", "cerrada_con_pendiente_declarado (E(AP) ICRP-116)"),
                          ("F77", "cerrada_en_F85 (mallas v0.35)"), ("F78", "sustituida_por_exterior_directo_F85")):
            f.write(json.dumps({"fecha": out["fecha"], "fase_afectada": fase, "estado": est}, ensure_ascii=False) + "\n")
    with open(OUT / "registro_mc.jsonl", "a") as f:
        f.write(json.dumps({"fecha": out["fecha"], "fase": 85, "md5_mallas": mallas["md5"], "historias": R2["historias_totales"]},
                           ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    o = ejecutar()
    print("F70 fuera bloque:", o["F70"]["fuera_bloque_rel"], " TVL medio", o["F71"]["tvl_medio_cm"])
    for r in o["F71"]["tabla_muros"]:
        print(r)
    for k, m in o["F78"]["muros"].items():
        print(k, f"ext {m['exterior_cara_uSv_h']:.3g}±{m['exterior_cara_se']:.2g}  parche {m['exterior_parche_uSv_h']:.3g}  anual {m['dosis_anual_mSv'][1]:.3g} mSv  ok={m['cumple_publico'] and m['fuera_de_zona']}")
    print({k: f"{v['A_sat_Bq']:.3g}±{(v['se_Bq'] or 0):.2g}" for k, v in o["F72"].items()})
    print("F77:", {k: v for k, v in o["F77"].items() if k in ("voxeles", "finitas", "no_negativas", "voxeles_con_dosis", "voxeles_std_menor_50pct")})
