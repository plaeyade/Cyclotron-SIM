#!/usr/bin/env python3
"""F82 — Re-corrida de F66 (TPMC), F74 (analítico Chiggiato) y F65 (stripping en gas) con la cámara y el bombeo del
gemelo kiube v0.35 (marcas de afectación de F81). Las evidencias de F66/F74/F65 NO se tocan: esta fase escribe la suya
y deja la comparación V04 → v0.35 cuantificada.

Entradas (kiube v0.35, verificadas por md5 contra el manifiesto de F81):
  - gap efectivo de cámara 102,6 mm (validation/gap_efectivo_camara.json; V04: 132);
  - 4 bombas ODP en los agujeros de valle (r = 127 mm a r = 375 mm, T04F2): 2 difusoras bajo la máquina (az 0° y 180°,
    tapa inferior; renders DS) + 2 en las carcasas RF (az 90° y 270°, tapa superior; ubicación EST). V04: 4 bocas a 45°
    (sobre las colinas) en ambas tapas;
  - bombeo por bomba: brida (H₂O 1300 l/s, H₂ 2000 l/s; EST) en serie con el agujero A+C (Clausing 0,42, CALC);
  - cargas de gas: H₂O de paredes y juntas a 10 h 5,2e-4 mbar·l/s (CERN-CH con acero ×30 calibrado al bombeo medido de
    TUD03T) y H₂ de la fuente 2 sccm = 3,4e-2 mbar·l/s (EST; AMIT usa 3,5 sccm). V04: 2,3e-6 + 1e-3 mbar·l/s.
Dos gases → dos corridas TPMC (H₂O m = 18, H₂ m = 2,016) y suma de presiones. Stripping: σ ∝ 1/β² para todos los gases
(PLUM-CYR), así que la mezcla se reduce a una presión H₂-equivalente P_H2 + (σ_Nakhai/σ_H2)·P_H2O.
Contraste independiente: transmisión fuente→stripper medida por IBA 81 % (TUD03), de la que kiube atribuye 15 % al gas.
Reglas: simulacion_teorica, pq no_declarada_nunca, append-only, 18.10, memoria 3 (sin haz real).
Rutas: KIUBE_RAIZ, RPOS_EVIDENCIA.
"""
from __future__ import annotations
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import datetime, hashlib, json, math, os, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]
from app.domain.vacuum_tpmc import CamaraVacio, simular_mapa_presion, R_GAS        # noqa: E402
from app.domain.gas_stripping import (PerfilPresion, simular_perdida_gas, sigma_mezcla_cm2,   # noqa: E402
                                      sigma_molecula_cm2, COMPOSICION_RESIDUAL_NAKHAI)

KIUBE = pathlib.Path(os.environ.get("KIUBE_RAIZ", f"{_RPOS_RAIZ}/work/kiube_v035/kiube"))
EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072"))
OUT = EV / "fase82"
N_MOL = int(os.environ.get("F82_N_MOL", "600000"))
T_K = 293.15
SCCM = 0.016887                   # mbar·l/s por sccm a 0 °C/1013 mbar (definición)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def v_media_cm_s(m_g_mol):
    return math.sqrt(8 * R_GAS * T_K / (math.pi * m_g_mol / 1000.0)) * 100.0


def entradas():
    import yaml
    man = {a["ruta"]: a["md5"] for a in json.load(open(KIUBE / "out/rpos_v035/manifiesto.json"))["artefactos"]}
    usados = ["params.yaml", "validation/gap_efectivo_camara.json", "out/vacuum/report.json"]
    for r in usados:
        assert md5(KIUBE / r) == man[r], f"md5 de {r} no coincide con el manifiesto F81"
    P = yaml.safe_load(open(KIUBE / "params.yaml"))
    gap = json.load(open(KIUBE / "validation/gap_efectivo_camara.json"))
    rep = json.load(open(KIUBE / "out/vacuum/report.json"))["presion"]
    vac = P["vacuum"]
    return dict(
        md5={r: man[r] for r in usados},
        gap_mm=gap["gap_efectivo_mm"], r_cam=P["vacuum_chamber"]["inner_radius"]["v"],
        r_boca=P["yoke"]["pump_hole_radius"]["v"], r_pos=P["yoke"]["pump_hole_r_pos"]["v"],
        W=vac["port_transmission"]["v"], S_brida={"H2O": vac["pump_speed_H2O_l_s"]["v"], "H2": vac["pump_speed_l_s"]["v"]},
        Q_H2O=rep["Q_H2O_10h_mbar_l_s"], Q_H2=vac["h2_sccm_per_active_source"]["v"] * SCCM,
        sccm=vac["h2_sccm_per_active_source"]["v"], S_kiube={"H2O": rep["S_eff_H2O_l_s"], "H2": rep["S_eff_H2_l_s"]},
        P_kiube={"H2O_10h": rep["P_H2O_10h_mbar"], "H2_operacion": rep["P_H2_operacion_mbar"]},
        perdida_gas_kiube=vac["gas_loss_fraction"]["v"], P_base_medida=vac["P_base_measured_mbar"]["v"])


def bombas(E, gas, m):
    A = math.pi * (E["r_boca"] / 10) ** 2
    C_ap = A * v_media_cm_s(m) / 4 / 1000.0                          # l/s, apertura ideal
    S_p = 1.0 / (1.0 / E["S_brida"][gas] + 1.0 / (C_ap * E["W"]))   # brida en serie con el agujero A+C
    s = S_p / C_ap
    lst = []
    for az, lado, quien in ((0, -1, "difusora bajo la máquina"), (180, -1, "difusora bajo la máquina"),
                            (90, +1, "carcasa RF (EST)"), (270, +1, "carcasa RF (EST)")):
        lst.append(dict(x_mm=E["r_pos"] * math.cos(math.radians(az)), y_mm=E["r_pos"] * math.sin(math.radians(az)),
                        lado=lado, r_mm=E["r_boca"], s=s, az_deg=az, ubicacion=quien))
    return lst, dict(C_apertura_l_s=C_ap, S_bomba_l_s=S_p, s_captura=s, S_total_l_s=4 * S_p)


def ejecutar():
    OUT.mkdir(parents=True, exist_ok=True)
    E = entradas()
    res = {}
    for gas, m, Q in (("H2O", 18.015, E["Q_H2O"]), ("H2", 2.016, E["Q_H2"])):
        b, an = bombas(E, gas, m)
        area = 2 * math.pi * (E["r_cam"] / 10) ** 2 + 2 * math.pi * (E["r_cam"] / 10) * (E["gap_mm"] / 10)
        cam = CamaraVacio(r_pared_mm=E["r_cam"], gap_mm=E["gap_mm"], n_bombas=4, r_bomba_mm=E["r_boca"],
                          r_pos_bomba_mm=E["r_pos"], s_captura=an["s_captura"],
                          q_pared_mbar_l_s_cm2=(Q / area if gas == "H2O" else 0.0),
                          q_h2_fuentes_total_mbar_l_s=(Q if gas == "H2" else 0.0),
                          masa_gas_g_mol=m, temperatura_k=T_K,
                          bombas=tuple({k: v for k, v in x.items() if k in ("x_mm", "y_mm", "lado", "r_mm", "s")} for x in b))
        mapa = simular_mapa_presion(cam, n_moleculas=N_MOL, n_bins_r=20)
        mapa.pop("fuentes", None)
        mapa["camara"] = {k: v for k, v in mapa["camara"].items() if k != "bombas"}
        res[gas] = dict(mapa=mapa, bombas=b, analitico=an, Q_mbar_l_s=Q,
                        P_eq_analitica_mbar=Q / an["S_total_l_s"],
                        P_media_tpmc_mbar=sum(mapa["p_mbar"]) / len(mapa["p_mbar"]),
                        captura_tpmc=mapa["absorbidos"] / max(mapa["hits_bomba"], 1))
    r_mm = res["H2"]["mapa"]["r_mm"]
    P_tot = [a + b for a, b in zip(res["H2O"]["mapa"]["p_mbar"], res["H2"]["mapa"]["p_mbar"])]
    # F65 con el perfil v0.35: presión H2-equivalente (σ ∝ 1/β² para todos los gases)
    k_sig = sigma_mezcla_cm2(COMPOSICION_RESIDUAL_NAKHAI, 18.0) / sigma_molecula_cm2("H2", 18.0)
    P_eqH2 = [h + k_sig * w for h, w in zip(res["H2"]["mapa"]["p_mbar"], res["H2O"]["mapa"]["p_mbar"])]
    r_ext = 463.3                                                     # órbita de 18 MeV en la lámina (kiube v0.35)
    B_ef = math.sqrt(18.0 ** 2 + 2 * 18.0 * 939.294) / 299.792458 / (r_ext / 1000)
    ESC = dict(energia_final_mev=18.0, delta_e_kev_vuelta=4 * 40.0 * math.sin(math.radians(39.0)), b_media_t=B_ef,
               r_pared_mm=E["r_cam"])
    perfil = PerfilPresion.desde_mapa(r_mm, P_eqH2)
    stripping = {}
    for nombre, I in (("KIUBE180_180uA", 180.0), ("comparable_F65_200uA", 200.0)):
        r = simular_perdida_gas(corriente_ua=I, perfil=perfil, composicion={"H2": 1.0}, **ESC)
        stripping[nombre] = dict(fraccion_perdida=r.fraccion_perdida, pct=round(100 * r.fraccion_perdida, 4),
                                 potencia_perdida_w=round(r.potencia_perdida_w, 2),
                                 energia_media_h0_mev=round(r.parametros["energia_media_h0_mev"], 2),
                                 primera_vuelta_pct=round(100 * r.perdida_por_vuelta[0], 4))
    # sensibilidad: H₂ de la fuente 3,5 sccm (AMIT)
    f_amit = 3.5 / E["sccm"]
    perfil_a = PerfilPresion.desde_mapa(r_mm, [f_amit * h + k_sig * w for h, w in
                                                zip(res["H2"]["mapa"]["p_mbar"], res["H2O"]["mapa"]["p_mbar"])])
    ra = simular_perdida_gas(corriente_ua=180.0, perfil=perfil_a, composicion={"H2": 1.0}, **ESC)
    # referencias V04 (evidencias anteriores, si existen en este entorno)
    v04 = dict(S_eff_l_s=7559.4, P_eq_F74_mbar=1.33e-7, perdida_F66_pct=0.00123, gap_mm=132.0, Q_total_mbar_l_s=1.002e-3)
    gas_ = stripping["KIUBE180_180uA"]["fraccion_perdida"]
    out = {
        "fase": 82, "version": "3.8.0", "tipo": "simulacion_teorica", "pq": "no_declarada_nunca",
        "fecha": datetime.datetime.now().isoformat(timespec="seconds"),
        "titulo": "Re-corrida F66/F74/F65 con cámara y bombeo kiube v0.35",
        "entradas": {k: v for k, v in E.items()}, "n_moleculas_por_gas": N_MOL, "temperatura_k": T_K,
        "F66_tpmc": {g: {"P_mbar_por_r": res[g]["mapa"]["p_mbar"], "P_media_mbar": res[g]["P_media_tpmc_mbar"],
                         "P_centro_mbar": res[g]["mapa"]["p_mbar"][0], "P_pared_mbar": res[g]["mapa"]["p_mbar"][-1],
                         "absorbidos": res[g]["mapa"]["absorbidos"], "hits_bomba": res[g]["mapa"]["hits_bomba"],
                         "captura_media_tpmc": res[g]["captura_tpmc"], "bombas": res[g]["bombas"]} for g in res},
        "r_mm": r_mm, "P_total_mbar": P_tot,
        "F74_analitico": {g: {**res[g]["analitico"], "Q_mbar_l_s": res[g]["Q_mbar_l_s"],
                              "P_eq_Q_S_mbar": res[g]["P_eq_analitica_mbar"],
                              "tpmc_media_sobre_analitica": res[g]["P_media_tpmc_mbar"] / res[g]["P_eq_analitica_mbar"]}
                          for g in res},
        "consistencia_con_kiube": {
            "S_H2O_l_s": [res["H2O"]["analitico"]["S_total_l_s"], E["S_kiube"]["H2O"]],
            "S_H2_l_s": [res["H2"]["analitico"]["S_total_l_s"], E["S_kiube"]["H2"]],
            "P_H2O_10h_mbar": [res["H2O"]["P_media_tpmc_mbar"], E["P_kiube"]["H2O_10h"]],
            "P_H2_operacion_mbar": [res["H2"]["P_media_tpmc_mbar"], E["P_kiube"]["H2_operacion"]],
            "nota": "kiube v0.35 usa su TPMC 3D con la geometría CAD completa; F82 usa el motor cilíndrico de F66 con gap efectivo"},
        "F65_stripping": {**stripping, "escenario": ESC, "sigma_mezcla_sobre_sigma_H2": k_sig,
                          "sensibilidad_H2_3p5sccm_AMIT_pct": round(100 * ra.fraccion_perdida, 4)},
        "comparacion_V04": {
            "S_eff_H2O_l_s": [v04["S_eff_l_s"], res["H2O"]["analitico"]["S_total_l_s"]],
            "P_eq_H2O_mbar": [v04["P_eq_F74_mbar"], res["H2O"]["P_eq_analitica_mbar"]],
            "perdida_gas_pct": [v04["perdida_F66_pct"], stripping["comparable_F65_200uA"]["pct"]],
            "perdida_V04_corregida_unidades_pct": v04["perdida_F66_pct"] * 1000.0,
            "factor_por_unidades": 1000.0,
            "factor_por_geometria_bombeo_y_carga_v035": stripping["comparable_F65_200uA"]["pct"] / (v04["perdida_F66_pct"] * 1000.0),
            "P_H2O_tpmc_sobre_medida_TUD03T": res["H2O"]["P_media_tpmc_mbar"] / E["P_base_medida"]},
        "contraste_IBA": {
            "transmision_medida_TUD03": 0.81, "perdida_gas_atribuida_kiube": E["perdida_gas_kiube"],
            "perdida_gas_F82": gas_, "cociente_F82_sobre_kiube": gas_ / E["perdida_gas_kiube"],
            "P_base_medida_TUD03T_mbar": E["P_base_medida"], "P_H2O_F82_mbar": res["H2O"]["P_media_tpmc_mbar"],
            "nota": ("La carga de H₂ (2 sccm) es EST y no se calibró con el 81 %: el contraste es independiente en la "
                     "magnitud. σ PLUM-CYR y órbita 2D de F65 sin cambios.")},
        "correccion_unidades": ("gas_stripping.MBAR_A_PA corregido de 0,1 a 100 (1 mbar = 100 Pa): todas las pérdidas por "
                                "gas de F65/F66 eran 1000× bajas; registro en fase82/correcciones_historicas.jsonl"),
        "advertencias": [
            "Simulación teórica: sin haz real (memoria 3). No cierra 18.19; PQ no_declarada_nunca (14.9).",
            "Motor cilíndrico de F66 con gap efectivo: no resuelve colinas (24 mm) ni valles (200 mm) por separado.",
            "Ubicación de las 2 bombas en carcasas RF: EST (TUD03 confirma 4 ODP, no su posición).",
            "Las evidencias de F66/F74/F65 no se modifican; F67 (Molflow+) sigue afectada y pendiente.",
        ],
    }
    json.dump(out, open(OUT / "vacio_v035_f82.json", "w"), indent=1, ensure_ascii=False, default=float)
    with open(OUT / "registro_vacio.jsonl", "a") as f:
        f.write(json.dumps({"fecha": out["fecha"], "fase": 82, "n_mol": N_MOL,
                            "P_media_H2O": res["H2O"]["P_media_tpmc_mbar"], "P_media_H2": res["H2"]["P_media_tpmc_mbar"],
                            "perdida_gas_180uA": gas_}, ensure_ascii=False) + "\n")
    # corrección histórica de unidades (append-only, precedente F73): las evidencias de F65/F66 no se reescriben
    corr = [dict(fase="F65", magnitud="fraccion_perdida (todos los escenarios)", factor=1000.0,
                 valor_registrado="0,0006–0,0063 % en 1e-7–1e-6 mbar", valor_corregido="0,6–6,3 %",
                 causa="gas_stripping.MBAR_A_PA = 0,1 (1 mbar = 100 Pa): densidad n = P/kT 1000× baja"),
            dict(fase="F66", magnitud="pérdida H⁻ con mapa TPMC", factor=1000.0, valor_registrado="0,00123 %",
                 valor_corregido="≈1,23 %", causa="misma constante; el mapa de presión de F66 es correcto"),
            dict(fase="F66/F67/F74", magnitud="presiones", factor=1.0, valor_registrado="sin cambio",
                 valor_corregido="sin cambio", causa="vacuum_tpmc: dos factores de unidades se compensaban exactamente")]
    with open(OUT / "correcciones_historicas.jsonl", "a") as f:
        for c in corr:
            f.write(json.dumps({"fecha": out["fecha"], **c}, ensure_ascii=False) + "\n")
    # resolución de las marcas de F81 (append-only; la marca original no se toca)
    with open(OUT / "resolucion_marcas_f81.jsonl", "a") as f:
        for fase in ("F66", "F74", "F65"):
            f.write(json.dumps({"fecha": out["fecha"], "fase_afectada": fase, "estado": "re_corrida_en_F82",
                                "evidencia": "fase82/vacio_v035_f82.json"}, ensure_ascii=False) + "\n")
        f.write(json.dumps({"fecha": out["fecha"], "fase_afectada": "F67", "estado": "afectada_pendiente_re_corrida",
                            "motivo": "requiere Molflow+ oficial con la geometría v0.35"}, ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    o = ejecutar()
    for g in ("H2O", "H2"):
        t, a = o["F66_tpmc"][g], o["F74_analitico"][g]
        print(f"{g}: S = {a['S_total_l_s']:.0f} l/s (s = {a['s_captura']:.3f}), P_eq Q/S = {a['P_eq_Q_S_mbar']:.2e}; "
              f"TPMC centro {t['P_centro_mbar']:.2e} → pared {t['P_pared_mbar']:.2e}, media {t['P_media_mbar']:.2e} mbar")
    print("stripping:", {k: v for k, v in o["F65_stripping"].items() if k.startswith(("KIUBE", "comparable", "sens"))})
    print("V04 →", o["comparacion_V04"])
    print("IBA:", {k: o["contraste_IBA"][k] for k in ("perdida_gas_atribuida_kiube", "perdida_gas_F82", "cociente_F82_sobre_kiube")})
