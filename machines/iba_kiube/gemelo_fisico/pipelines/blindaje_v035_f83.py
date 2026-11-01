#!/usr/bin/env python3
"""F83 — Re-corrida de F70 (MC máquina + búnker BIND) con la geometría kiube v0.35 — marca de afectación de F81.

Motor: Geant4 11.4 (QGSP_BIC_HP, G4NDL) en lugar de OpenMC 0.16: en este entorno no hay biblioteca ENDF/B-VIII.0 en
HDF5 y la red solo permite registros de paquetes/GitHub. La evidencia de F70 (OpenMC) no se toca; esta fase es un
cálculo independiente con otro código (lo que además sirve de contraste cruzado de F70).
Geometría: máquina kiube v0.35 teselada (GDML del LOD neutrónico; blancos F18 T1/T5 en la cara del yugo, r ≈ 969 mm,
con puertas de HDPE de 18 cm) en el bloque BIND igual al de F70 (sala 5,9 × 4,7 × 2,5 m, muros 2,11 m, losa 0,55 m,
techo 1,45 m). Fuente 18O(p,n) BON21 (4,93e-3 n/p), 2 × 90 µA (KIUBE 180, PROP).
Reducción de varianza: cascarones laterales de importancia cada 10 cm (factor 2,15) a través de los muros → cuentas
directas fuera del muro (F70: 0 cuentas; F78 las extrapoló).
Entradas (repositorio kiube): out/mc_f70v035/<caso>/resultado_mapa.json y mapa_h10.npz (md5 registrado).
Reglas: simulacion_teorica, pq no_declarada_nunca, append-only, 18.10, memoria 3.
Rutas: KIUBE_RAIZ, RPOS_EVIDENCIA, F83_CASO (por defecto out/mc_f70v035/h2500_t1450).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", f"{_RPOS_RAIZ}")  # raíz de evidencia/datos externos
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", (_RPOS_RAIZ + ""))  # raíz de evidencia/datos externos
import datetime, hashlib, json, math, os, pathlib

KIUBE = pathlib.Path(os.environ.get("KIUBE_RAIZ", (_RPOS_RAIZ + "/work/kiube_v035/kiube")))
CASO = KIUBE / os.environ.get("F83_CASO", "out/mc_f70v035/h2500_t1450")
EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", (_RPOS_RAIZ + "/work/evidencia_v072")))
OUT = EV / "fase83"

# valores de referencia de F70/F72/F78 citados de los mensajes de commit/CHANGELOG del repositorio (la evidencia JSON de
# esas fases vive en evidencia_v072, que no está en este entorno: se citan, no se recalculan)
REF = {
    "F70_fuera_bloque_rel_max": 0.05,                 # test_mc_f70: 0 < dose_integral_fuera_bloque_rel < 0.05
    "F70_exterior_cuentas": 0,                        # F76/F78: exterior en 0 cuentas
    "F78_mu_cm_inv": 0.0142,                          # commit F78 v3.4.0
    "F78_T_H10": 4.13e-2,                             # commit F78: transmisión estimada interior 287 → exterior 512 cm
    "F72_A_sat_Bq": {"Co60": 6.4e9, "Eu152": 6.4e9, "Ag110m": 23e9, "Cu64": 861e9},   # commit F73 (unidades corregidas)
}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def ejecutar():
    OUT.mkdir(parents=True, exist_ok=True)
    R = json.load(open(CASO / "resultado_mapa.json"))
    caras = R["caras_plano_medio"]
    # transmisión directa del muro: exterior (0–30 cm fuera, media de la cara) / interior (capa de sala junto al muro)
    trans, mus = {}, {}
    for k, c in caras.items():
        ext, se = c["exterior_0_30cm_cara_uSv_h"]
        trans[k] = dict(exterior_uSv_h=ext, exterior_se=se, interior_uSv_h=c["interior_cara_uSv_h"],
                        T_directa=(ext / c["interior_cara_uSv_h"]) if c["interior_cara_uSv_h"] else None,
                        mu_cm_inv=c["mu_cara_cm_inv"], tvl_cm=c["tvl_cara_cm"])
        if c["mu_cara_cm_inv"]:
            mus[k] = c["mu_cara_cm_inv"]
    mu_med = sum(mus.values()) / len(mus) if mus else None
    T_225_v035 = math.exp(-mu_med * 225) if mu_med else None
    act = R["activacion_saturacion"]
    comp_act = {iso: dict(F83_Bq=act[iso]["A_sat_Bq"], F83_se_Bq=act[iso]["sigma_Bq"], F72_Bq=REF["F72_A_sat_Bq"][iso],
                          cociente_F83_F72=act[iso]["A_sat_Bq"] / REF["F72_A_sat_Bq"][iso], fuente_xs=act[iso]["fuente_xs"])
                for iso in act}
    out = {
        "fase": 83, "version": "3.9.0", "tipo": "simulacion_teorica", "pq": "no_declarada_nunca",
        "fecha": datetime.datetime.now().isoformat(timespec="seconds"),
        "titulo": "Re-corrida F70 con geometría kiube v0.35 (Geant4, búnker BIND, blancos en cara del yugo con puertas HDPE)",
        "motor": {"codigo": "Geant4 11.4.3", "fisica": "QGSP_BIC_HP (G4NDL)", "motivo_cambio": (
            "sin ENDF/B-VIII.0 HDF5 para OpenMC en este entorno; F70 (OpenMC) queda como evidencia y esta fase es un "
            "cálculo independiente con otro código")},
        "entradas": {"caso": str(CASO.relative_to(KIUBE)), "md5_resultado": md5(CASO / "resultado_mapa.json"),
                     "md5_mapa": md5(CASO / "mapa_h10.npz"), "corridas": R["corridas"], "historias": R["historias_totales"],
                     "fuente": R["fuente"], "bloque": R["bloque_hormigon"]},
        "resultados": {
            "integral_dosis_fuera_bloque_rel": R["integral_dosis_fuera_bloque_rel"],
            "puntos_1m_T1_uSv_h": R["srg_tabla8_uSv_h"],
            "muros": trans, "mu_medio_cm_inv": mu_med, "tvl_medio_cm": math.log(10) / mu_med if mu_med else None,
            "techo_exterior_medio_uSv_h": R["techo_exterior_medio_uSv_h"],
            "activacion": comp_act},
        "comparacion": {
            "F70_fuera_bloque_rel": {"F70_criterio": f"< {REF['F70_fuera_bloque_rel_max']}", "F83": R["integral_dosis_fuera_bloque_rel"],
                                     "cumple": R["integral_dosis_fuera_bloque_rel"] < REF["F70_fuera_bloque_rel_max"]},
            "exterior_directo": {"F70": "0 cuentas (F76/F78)", "F83": {k: v["exterior_uSv_h"] for k, v in trans.items()}},
            "F78_extrapolacion": {
                "mu_F78_cm_inv": REF["F78_mu_cm_inv"], "mu_F83_cm_inv": mu_med,
                "T_225cm_F78": REF["F78_T_H10"], "T_225cm_F83": T_225_v035,
                "cociente_T_F78_sobre_F83": (REF["F78_T_H10"] / T_225_v035) if T_225_v035 else None,
                "veredicto": "discrepancia: la μ de F78 (TVL 162 cm) no es física para hormigón de 2,35 g/cm³; el SRG "
                             "(Tabla 9: 200 → 220 cm, ×4,6) da TVL ≈ 30 cm, y F83 23–29 cm. La dosis exterior de F78 "
                             "está sobreestimada en muchos órdenes de magnitud (conservadora, pero no representativa)."},
            "F72_activacion": "mismo orden en Co60/Eu152/Cu64; Ag110m difiere (I0 de Ag109→110m EST en F83; F72 sin desglose aquí). "
                              "Normalización de corriente de F72 no disponible en este entorno: los cocientes son indicativos.",
            "SRG_tabla8": "los puntos a 1 m incluyen la puerta de HDPE (configuración del proyecto); la Tabla 8 del SRG es sin "
                          "puertas: no son comparables directamente (el benchmark sin puertas está en kiube v0.27, 0,51–1,22)."},
        "advertencias": [
            "Simulación teórica: sin haz real (memoria 3). No cierra 18.19; PQ no_declarada_nunca (14.9).",
            "Código distinto de F70 (Geant4 frente a OpenMC): las diferencias incluyen geometría (v0.8 → v0.35) y código.",
            "Exterior: estadística directa con error estándar entre corridas (cascarones de importancia); profundidades > 1,5 m con error > 50 %.",
            "Artefacto de varianza declarado: entre 1,4 y 2,0 m de profundidad el perfil medio forma una meseta con barras de error grandes "
            "(neutrones que cruzan el techo sin sesgo y entran en cascarones exteriores del muro con importancia alta). La esperanza es "
            "insesgada pero la varianza es enorme; μ se ajusta en 40–130 cm y el exterior (0–30 cm fuera) es coherente con extrapolar la "
            "exponencial. Corrección pendiente: limitar los cascarones a la altura de la sala.",
            "Activación: convención de Høgdahl (σ0·φ_th + I0·φ_epi/u); σ de Eu151→152g e I0 de Ag109→110m son EST de orden de tabla.",
            "Ventana de techo 1,45 m y losa 0,55 m: iguales a F70 por comparabilidad (techo real pendiente de los cortes; SRG mín. 1,7 m).",
        ],
    }
    json.dump(out, open(OUT / "blindaje_v035_f83.json", "w"), indent=1, ensure_ascii=False, default=float)
    with open(OUT / "registro_mc.jsonl", "a") as f:
        f.write(json.dumps({"fecha": out["fecha"], "fase": 83, "md5_mapa": out["entradas"]["md5_mapa"],
                            "historias": R["historias_totales"], "fuera_bloque_rel": R["integral_dosis_fuera_bloque_rel"]},
                           ensure_ascii=False) + "\n")
    with open(OUT / "resolucion_marcas_f81.jsonl", "a") as f:
        for fase, est, nota in (("F70", "re_corrida_en_F83", "Geant4, geometría v0.35"),
                                ("F72", "re_corrida_parcial_en_F83", "activación de saturación por región"),
                                ("F78", "contradicha_por_F83", "μ 0,0142 cm⁻¹ frente a ~0,09 cm⁻¹ directo"),
                                ("F71/F76/F77", "afectada_pendiente_re_corrida", "comparan con F70 v0.8")):
            f.write(json.dumps({"fecha": out["fecha"], "fase_afectada": fase, "estado": est, "nota": nota}, ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    o = ejecutar()
    r = o["resultados"]
    print("fuera bloque:", r["integral_dosis_fuera_bloque_rel"], " μ medio", r["mu_medio_cm_inv"], " TVL", r["tvl_medio_cm"])
    for k, v in r["muros"].items():
        print(k, {a: (f"{b:.3g}" if isinstance(b, float) else b) for a, b in v.items()})
    print("F78:", o["comparacion"]["F78_extrapolacion"]["cociente_T_F78_sobre_F83"])
    print({k: round(v["cociente_F83_F72"], 3) for k, v in r["activacion"].items()})
