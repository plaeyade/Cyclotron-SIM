#!/usr/bin/env python3
"""F84 — Re-corrida de F67 (validación cruzada con Molflow+ oficial) con la cámara y el bombeo de kiube v0.35 (F82).

Cierra la marca de afectación de F81 sobre F67. Mismo método que F67: tools/kiube_camara_molflow.py genera la cámara
(ahora con la configuración de F82: gap 102,6 mm, 2 bocas en la tapa inferior a 0°/180° y 2 en la superior a 90°/270°,
r = 127 mm a 375 mm, captura de F82 por gas) y molflowCLI 2.11.1 (binario oficial CERN, paquete Debian del usuario)
calcula la presión por faceta. Dos gases, como F82: H₂O de paredes (5,2e-4 mbar·l/s) y H₂ de la fuente (2 sccm).
Se compara con el mapa TPMC de F82 y se recalcula la pérdida por gas con el mapa de Molflow.
Validación analítica previa (ambos motores): cilindro r = 50 cm, h = 100 cm, una boca r = 5 cm con s = 1 → P = Q/S
uniforme. Encontró dos errores que se corrigen aquí (registro append-only, evidencias anteriores intactas):
  - vacuum_tpmc (F66/F82): área de una sola tapa para los impactos de las dos → presiones ×2;
  - tools/kiube_camara_molflow.py (F67): la tapa terminaba en el último anillo de 3 cm → corona abierta entre tapa y
    pared, el gas escapaba → presiones de Molflow bajas (×16 en el caso analítico);
  - tools/kiube_camara_molflow.py (F67): la fuente de H₂ (disco de 2 caras) inyectaba la mitad de su caudal.
Con ambas correcciones los dos motores reproducen Q/S al 0,1–0,2 %. F82 se recalcula aquí con el motor corregido.
Reglas: simulacion_teorica, pq no_declarada_nunca, append-only, 18.10, memoria 3.
Rutas: RPOS_EVIDENCIA, MOLFLOW_CLI (binario), F84_NDES (desorciones por gas; por defecto 2e6).
"""
from __future__ import annotations
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import csv, datetime, hashlib, json, math, os, pathlib, subprocess, sys
import xml.etree.ElementTree as ET
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS); sys.path.insert(0, str(ROOT / "tools"))
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]
from app.domain.gas_stripping import PerfilPresion, simular_perdida_gas            # noqa: E402
import kiube_camara_molflow as KM                                                  # noqa: E402

EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072"))
OUT = EV / "fase84"
CLI = os.environ.get("MOLFLOW_CLI", "molflowCLI")
NDES = float(os.environ.get("F84_NDES", "2e6"))
KB = 1.380649e-23


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def correr(xml, tag):
    res = OUT / f"res_{tag}"
    subprocess.run([CLI, "-f", str(xml), "-d", f"{NDES:g}", "-j", "2", "--writeFacetQuantities", "--outputPath", str(res),
                    "-o", f"{tag}_out.zip", "--noProgress", "-s", "0", "-a", "0"], check=True,
                   stdout=open(OUT / f"{tag}.log", "w"), stderr=subprocess.STDOUT)
    return res / "facet_physics.csv"


def leer(csv_path, xml_path, gap_mm, r_mm=548.0, n_bins=20, temp=293.15):
    root = ET.parse(xml_path).getroot()
    V = np.array([[float(v.get(k)) for k in "xyz"] for v in root.find("Geometry/Vertices")])
    rows = list(csv.DictReader(open(csv_path)))
    P = np.array([float(r["Pressure [mbar]"]) for r in rows]); Zt = np.array([float(r["Imping.rate"]) for r in rows])
    rad, zz, area, st = [], [], [], []
    for f in root.find("Geometry/Facets"):
        idx = [int(i.get("vertex")) for i in f.find("Indices")]
        vs = V[idx]; c = vs.mean(0)
        rad.append(math.hypot(c[0], c[1]) * 10); zz.append(c[2] * 10)
        area.append(sum(0.5 * np.linalg.norm(np.cross(vs[k] - vs[0], vs[k + 1] - vs[0])) for k in range(1, len(vs) - 1)))
        st.append(float(f.find("Sticking").get("constValue")))
    rad, zz, area, st = map(np.array, (rad, zz, area, st))
    tapa = (np.abs(zz) > 0.4 * gap_mm) & (st < 0.01)                   # facetas de tapa (sin las bocas)
    b = np.linspace(0, r_mm, n_bins + 1); rc, pc = [], []
    for i in range(n_bins):
        m = tapa & (rad >= b[i]) & (rad < b[i + 1])
        if m.any():
            rc.append(float((b[i] + b[i + 1]) / 2)); pc.append(float(np.average(P[m], weights=area[m])))
    bomba = st >= 0.01
    q_abs = float((Zt[bomba] * area[bomba] * 1e-4 * st[bomba]).sum() * KB * temp * 10.0)   # mbar·l/s
    return rc, pc, q_abs


def ejecutar():
    OUT.mkdir(parents=True, exist_ok=True)
    # --- validación analítica de los dos motores ---
    import importlib.util
    from app.domain.vacuum_tpmc import CamaraVacio, simular_mapa_presion
    Ra, Ha, rpa, Qa = 50.0, 100.0, 5.0, 1e-3
    area_a = 2 * math.pi * Ra * Ra + 2 * math.pi * Ra * Ha
    va = math.sqrt(8 * 8.314462618 * 293.15 / (math.pi * 0.018015)) * 100
    Sa = math.pi * rpa ** 2 * va / 4 / 1000
    xml_a = OUT / "analitico.xml"; KM.TEMP = 293.15
    KM.generar_xml(str(xml_a), dict(R_cm=Ra, gap_cm=Ha, temp=293.15, masa=18.015, bombas=[(25.0, 0.0, 1, rpa, 1.0)],
                                    q_pared=Qa / area_a, q_fuente=1e-30))
    Pm = np.array([float(r["Pressure [mbar]"]) for r in csv.DictReader(open(correr(xml_a, "analitico")))])
    ma = simular_mapa_presion(CamaraVacio(r_pared_mm=Ra * 10, gap_mm=Ha * 10, n_bombas=1, r_bomba_mm=rpa * 10,
                                          r_pos_bomba_mm=250, s_captura=1.0, q_pared_mbar_l_s_cm2=Qa / area_a,
                                          q_h2_fuentes_total_mbar_l_s=0.0, masa_gas_g_mol=18.015, temperatura_k=293.15,
                                          bombas=({"x_mm": 250.0, "y_mm": 0.0, "lado": 1, "r_mm": rpa * 10, "s": 1.0},)),
                              n_moleculas=200_000, n_bins_r=10)
    analitico = dict(geometria="cilindro r 50 cm, h 100 cm, 1 boca r 5 cm s 1, desgasificado uniforme 1e-3 mbar·l/s",
                     P_QS_mbar=Qa / Sa, P_molflow_mbar=float(Pm.mean()), P_tpmc_mbar=float(np.mean(ma["p_mbar"])),
                     molflow_sobre_QS=float(Pm.mean()) / (Qa / Sa), tpmc_sobre_QS=float(np.mean(ma["p_mbar"])) / (Qa / Sa),
                     antes_de_corregir=dict(molflow_sobre_QS=0.0635, tpmc_sobre_QS=2.003))
    # --- F82 recalculada con el motor corregido (en fase84/; la evidencia de F82 no se toca) ---
    spec = importlib.util.spec_from_file_location("f82", ROOT / "pipelines/vacio_v035_f82.py")
    F82 = importlib.util.module_from_spec(spec); spec.loader.exec_module(F82)
    F82.OUT = OUT / "f82_corregida"
    f82 = F82.ejecutar()
    f82_orig = json.load(open(EV / "fase82/vacio_v035_f82.json"))
    E = f82["entradas"]; gap = E["gap_mm"]
    res = {}
    # control: cámara V04 sin configuración (debe reproducir F67)
    xml_v04 = OUT / "control_v04.xml"; KM.TEMP = 293.15; KM.generar_xml(str(xml_v04))
    rc, pc, qa = leer(correr(xml_v04, "control_v04"), xml_v04, 132.0)
    mv = simular_mapa_presion(CamaraVacio(temperatura_k=293.15), n_moleculas=300_000, n_bins_r=20)   # V04 corregida
    pv = np.interp(rc, mv["r_mm"], mv["p_mbar"])
    control_tpmc = dict(P_centro_tpmc_mbar=float(mv["p_mbar"][0]), ratio_molflow_tpmc_medio=float(np.mean(np.array(pc) / pv)),
                        ratio_min=float(np.min(np.array(pc) / pv)), ratio_max=float(np.max(np.array(pc) / pv)))
    control = dict(**control_tpmc, P_centro_molflow_estanco_mbar=pc[0], P_centro_F67_registrado_mbar=1.25e-7, cociente=pc[0] / 1.25e-7,
                   nota="F67 tenía la corona abierta; F66 ponía bombas en ambas tapas y F67 solo en la superior (configuraciones distintas)")
    for gas, m, Q in (("H2O", 18.015, E["Q_H2O"]), ("H2", 2.016, E["Q_H2"])):
        b = f82["F66_tpmc"][gas]["bombas"]
        cfg = dict(R_cm=E["r_cam"] / 10, gap_cm=gap / 10, temp=f82["temperatura_k"], masa=m,
                   bombas=[(x["x_mm"] / 10, x["y_mm"] / 10, x["lado"], x["r_mm"] / 10, x["s"]) for x in b])
        area_cm2 = 2 * math.pi * (E["r_cam"] / 10) ** 2 + 2 * math.pi * (E["r_cam"] / 10) * (gap / 10)
        cfg["q_pared"] = Q / area_cm2 if gas == "H2O" else 0.0
        cfg["q_fuente"] = Q if gas == "H2" else 0.0
        if gas == "H2O":
            cfg["q_fuente"] = 1e-30                                       # disco de la fuente sin gas
        xml = OUT / f"camara_v035_{gas}.xml"
        info = KM.generar_xml(str(xml), cfg)
        rc, pc, qa = leer(correr(xml, f"v035_{gas}"), xml, gap, E["r_cam"], temp=f82["temperatura_k"])
        r82, p82 = f82["r_mm"], f82["F66_tpmc"][gas]["P_mbar_por_r"]
        p82i = np.interp(rc, r82, p82)
        ratio = np.array(pc) / p82i
        res[gas] = dict(xml=info, md5_xml=md5(xml), r_mm=rc, P_molflow_mbar=pc, P_tpmc_F82_mbar=p82i.tolist(),
                        ratio_medio=float(ratio.mean()), ratio_min=float(ratio.min()), ratio_max=float(ratio.max()),
                        P_media_molflow_mbar=float(np.mean(pc)), P_media_F82_mbar=float(np.mean(p82i)),
                        Q_inyectado_mbar_l_s=Q, Q_absorbido_mbar_l_s=qa, fraccion_absorbida=qa / Q if Q else None)
    # pérdida por gas con el mapa de Molflow (mismo escenario y presión H2-equivalente que F82)
    k = f82["F65_stripping"]["sigma_mezcla_sobre_sigma_H2"]
    r = res["H2"]["r_mm"]
    Peq = [h + k * w for h, w in zip(res["H2"]["P_molflow_mbar"], np.interp(r, res["H2O"]["r_mm"], res["H2O"]["P_molflow_mbar"]))]
    esc = f82["F65_stripping"]["escenario"]
    st = simular_perdida_gas(corriente_ua=180.0, perfil=PerfilPresion.desde_mapa(r, Peq), composicion={"H2": 1.0}, **esc)
    f82_loss = f82["F65_stripping"]["KIUBE180_180uA"]["fraccion_perdida"]
    out = {
        "fase": 84, "version": "3.10.0", "tipo": "simulacion_teorica", "pq": "no_declarada_nunca",
        "fecha": datetime.datetime.now().isoformat(timespec="seconds"),
        "titulo": "Re-corrida F67: Molflow+ 2.11.1 oficial con la cámara y el bombeo de kiube v0.35 (F82)",
        "motor": {"binario": "molflowCLI 2.11.1 (CERN, paquete Debian del usuario)", "md5_binario": md5(CLI) if os.path.exists(CLI) else None,
                  "desorciones_por_gas": NDES},
        "validacion_analitica": analitico,
        "control_v04": control,
        "f82_corregida": {"P_media_H2O_mbar": f82["F66_tpmc"]["H2O"]["P_media_mbar"], "P_media_H2_mbar": f82["F66_tpmc"]["H2"]["P_media_mbar"],
                          "P_media_H2O_F82_original_mbar": f82_orig["F66_tpmc"]["H2O"]["P_media_mbar"],
                          "perdida_gas_180uA": f82["F65_stripping"]["KIUBE180_180uA"]["fraccion_perdida"],
                          "perdida_gas_180uA_F82_original": f82_orig["F65_stripping"]["KIUBE180_180uA"]["fraccion_perdida"],
                          "perdida_gas_3p5sccm_pct": f82["F65_stripping"]["sensibilidad_H2_3p5sccm_AMIT_pct"],
                          "P_H2O_sobre_medida_TUD03T": f82["comparacion_V04"]["P_H2O_tpmc_sobre_medida_TUD03T"],
                          "evidencia": "fase84/f82_corregida/vacio_v035_f82.json"},
        "gases": res,
        "stripping": {"perdida_molflow": st.fraccion_perdida, "perdida_F82_tpmc": f82_loss,
                      "cociente_molflow_F82": st.fraccion_perdida / f82_loss,
                      "escenario": esc, "corriente_uA": 180.0},
        "observaciones": [
            "F67 calculó la pérdida con E = 17,4 MeV, ΔE = 80 keV/vuelta y B = 1,77 T, distinto del escenario de F65/F66 "
            "(18 MeV, 98,5 keV, 1,36 T): sus pérdidas no eran comparables entre fases. F84 usa el escenario de F82 (v0.35).",
            "Las pérdidas de F67 heredan además el factor 1000 de MBAR_A_PA (corregido en F82).",
        ],
        "advertencias": [
            "Simulación teórica: sin haz real (memoria 3). No cierra 18.19; PQ no_declarada_nunca (14.9).",
            "Misma cámara cilíndrica de gap efectivo que F82 (validación del motor, no de la geometría 3D completa).",
            "Ubicación de las 2 bombas en carcasas RF: EST.",
        ],
    }
    json.dump(out, open(OUT / "molflow_v035_f84.json", "w"), indent=1, ensure_ascii=False, default=float)
    with open(OUT / "registro_vacio.jsonl", "a") as f:
        f.write(json.dumps({"fecha": out["fecha"], "fase": 84, "ratio_H2O": res["H2O"]["ratio_medio"],
                            "ratio_H2": res["H2"]["ratio_medio"], "perdida_molflow": st.fraccion_perdida}, ensure_ascii=False) + "\n")
    with open(OUT / "correcciones_historicas.jsonl", "a") as f:
        for c in (dict(fase="F66/F82", magnitud="presiones TPMC", factor=0.5, causa="área de una sola tapa (vacuum_tpmc)"),
                  dict(fase="F65/F66/F82", magnitud="pérdida por gas con mapa TPMC", factor="≈0,5", causa="hereda la presión ×2"),
                  dict(fase="F67", magnitud="presiones Molflow", factor=control["cociente"],
                       causa="corona abierta tapa–pared en tools/kiube_camara_molflow.py (fuga de gas)"),
                  dict(fase="F67", magnitud="caudal de la fuente de H₂ en Molflow", factor=2.0,
                       causa="disco de 2 caras con desgasificación dividida entre 2 (Molflow la toma como total)"),
                  dict(fase="F66", magnitud="bombeo V04 del TPMC", factor="8 → 4 bocas",
                       causa="bocas activas en ambas tapas frente a las 4 de la tapa superior declaradas en F67/F74"),
                  dict(fase="F74", magnitud="P = Q/S analítica", factor=1.0, causa="no afectada (no usa los mapas)")):
            f.write(json.dumps({"fecha": out["fecha"], **c}, ensure_ascii=False, default=float) + "\n")
    with open(OUT / "resolucion_marcas_f81.jsonl", "a") as f:
        f.write(json.dumps({"fecha": out["fecha"], "fase_afectada": "F67", "estado": "re_corrida_en_F84",
                            "evidencia": "fase84/molflow_v035_f84.json"}, ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    o = ejecutar()
    print("analítico:", {k: v for k, v in o["validacion_analitica"].items() if k != "geometria"})
    print("control V04:", o["control_v04"])
    print("F82 corregida:", o["f82_corregida"])
    for g in ("H2O", "H2"):
        v = o["gases"][g]
        print(g, {k: (f"{x:.3g}" if isinstance(x, float) else x) for k, x in v.items() if k in
                  ("ratio_medio", "ratio_min", "ratio_max", "P_media_molflow_mbar", "P_media_F82_mbar", "fraccion_absorbida")})
    print("stripping:", {k: (f"{x:.4g}" if isinstance(x, float) else x) for k, x in o["stripping"].items() if k != "escenario"})
