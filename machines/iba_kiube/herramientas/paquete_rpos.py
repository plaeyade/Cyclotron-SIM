"""Paquete de ingesta para Radiopharma-OS (v0.35): artefactos multiformato con md5 y metadatos, verificación cruzada
contra los supuestos KIUBE-V04 / inventario F75 / SRG / KIUBE-OM, registro de re-etiquetado (append-only), impacto en
cadena y telemetría simulada declarada.

Reglas (las del repositorio RPOS):
  - nunca se sobrescribe un supuesto: el re-etiquetado es un registro append-only con el valor y la etiqueta anteriores;
  - re-etiquetar a `fuente_modelo_oficial_IBA (vía usuario)` exige DOS condiciones: que el valor coincida dentro de la
    tolerancia declarada Y que la cota del modelo tenga fuente IBA (documento oficial o publicación de autores IBA).
    Una coincidencia con una cota EST/CALC del modelo no confirma nada: queda `supuesto_declarado` con la nota;
  - toda diferencia se reporta cuantificada (Δ absoluto y relativo);
  - la telemetría del gemelo es `simulado_declarado`: no es haz real (memoria 3) y no puede entrar como medida al ajuste
    del Calibration Layer.
Uso: python herramientas/paquete_rpos.py [--out out/rpos_v035]
"""
from __future__ import annotations
import argparse, csv, datetime, hashlib, json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
VERSION = "v0.35"
ETQ_IBA = "fuente_modelo_oficial_IBA (vía usuario)"

# fuentes del modelo que son documentos IBA o publicaciones de autores IBA (docs/FUENTES.md)
FUENTES_IBA = {
    "SRG": "IBA Cyclone KIUBE Site Requirement Guide MID 55672 Rev. G (2023) — documento oficial",
    "PROP": "propuesta comercial IBA del proyecto — documento oficial",
    "DS": "ficha técnica IBA Cyclone KIUBE V6 — documento oficial",
    "DSR": "renders de la ficha técnica IBA — documento oficial (medida en imagen)",
    "T04F2": "TUP04 (Kleeven et al., IBA) Fig. 2, medida a escala ±3 mm — publicación IBA",
    "TUP04": "TUP04 (IBA) — publicación IBA", "TUD03": "TUD03 (IBA) — publicación IBA",
    "TUD03T": "TUD03, presentación (IBA) — publicación IBA", "KEY22": "Nuttens et al. (IBA), Cyclotrons2022 — publicación IBA",
    "VE24": "do Carmo et al. (IBA), Instruments 8 (2024) 38 — publicación IBA",
    "IPAC26": "Herrod et al. (IBA), IPAC2026 MOP7178 — publicación IBA", "PAT16": "patentes IBA 2016 — documento IBA",
}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def P_():
    import yaml
    return yaml.safe_load(open(ROOT / "params.yaml"))


def cota(raw, ruta):
    d = raw
    for k in ruta.split("."):
        d = d[k]
    return d


# ------------------------------------------------------------------------------------------------ 1. artefactos
REGLAS = [  # (patrón glob relativo, categoría, formato, unidades, generador)
    ("out/neutronics/**/*.step", "geometria", "STEP AP214", "mm", "build.py --lod neutronics"),
    ("out/fem_struct/**/*.step", "geometria", "STEP AP214", "mm", "build.py --lod fem_struct"),
    ("out/fem_mag/**/*.step", "geometria", "STEP AP214", "mm", "build.py --lod fem_mag"),
    ("out/print/*.step", "geometria", "STEP AP214", "mm", "build.py --lod print"),
    ("out/mc/*.h5m", "geometria", "DAGMC h5m (MOAB)", "cm", "build.py --lod mc"),
    ("out/mc/*.gdml", "geometria", "GDML (Geant4)", "mm", "build.py --lod mc"),
    ("out/mc/volumes.json", "metadatos_geometria", "JSON", "cm3", "build.py --lod mc"),
    ("out/print/stl/*.stl", "geometria", "STL binario", "mm (escala 1:1)", "build.py --lod print"),
    ("out/stl_hd/**/*.stl", "geometria", "STL binario", "mm (escala 1:10)", "build.py --lod stl_hd --escala 0.1"),
    ("out/stl_hd/**/*.3mf", "geometria", "3MF con color", "mm (escala 1:10)", "build.py --lod stl_hd --escala 0.1"),
    ("out/stl_hd/**/*.glb", "geometria", "glTF binario con color", "m (escala 1:10)", "build.py --lod stl_hd --escala 0.1"),
    ("out/vacuum/molflow/*", "geometria", "Molflow+ XML", "cm", "build.py --lod vacuum"),
    ("out/fem_mag/ngsolve/*.npz", "mapa_campo", "NumPy npz (r_m, phi_deg, Bz_T)", "m, °, T", "fem/magnet_ngsolve.py"),
    ("out/fem_mag/ngsolve/*.json", "curva", "JSON", "SI", "fem/magnet_ngsolve.py"),
    ("out/fem_rf/resultado_nominal.json", "curva", "JSON (modos RF)", "MHz, kW, J", "build.py --lod fem_rf"),
    ("out/orbitas/equilibrio_v34.json", "curva", "JSON (órbitas de equilibrio)", "MeV, m, MHz", "fem/orbitas.py"),
    ("out/vacuum/tpmc_pumpdown.json", "curva", "JSON (bombeo TPMC)", "h, mbar", "build.py --lod vacuum"),
    ("out/vacuum/report.json", "curva", "JSON", "mbar, l/s", "build.py --lod vacuum"),
    ("out/fuente_blanco/acel_p0_final.json", "curva", "JSON (aceleración 3D)", "MeV, m, °", "fem/fuente_a_blanco.py acel"),
    ("out/fuente_blanco/ext_p0_final.json", "curva", "JSON (extracción y mancha)", "mm, mrad, MeV", "fem/fuente_a_blanco.py ext"),
    ("validation/isocronia.json", "curva", "JSON", "%, MeV", "validation/figura_isocronia.py"),
    ("validation/fuente_blanco.json", "curva", "JSON", "%", "validation/figura_fuente_blanco.py"),
    ("validation/gap_efectivo_camara.json", "verificacion", "JSON", "m3, mm", "herramientas/gap_efectivo.py"),
    ("validation/validacion_srg.json", "verificacion", "JSON", "mm, kg", "validation/validar_srg.py"),
    ("validation/benchmark_srg.json", "verificacion", "JSON", "µSv/h", "validation/benchmark_srg.py"),
    ("params.yaml", "parametros", "YAML (cota + fuente + nota)", "mm, °, kA·v…", "fuente única del modelo"),
    ("materials.yaml", "parametros", "YAML (composición, densidad, curvas B-H)", "g/cm3, T, A/m", "fuente única del modelo"),
]


def artefactos():
    out, vistos = [], set()
    for pat, cat, fmt, uni, gen in REGLAS:
        for p in sorted(ROOT.glob(pat)):
            if not p.is_file() or p in vistos:
                continue
            vistos.add(p)
            st = p.stat()
            out.append(dict(ruta=str(p.relative_to(ROOT)), md5=md5(p), bytes=st.st_size, formato=fmt, categoria=cat,
                            unidades=uni, generador=gen, modelo="kiube", version_modelo=VERSION,
                            fecha_modificacion=datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
                            tipo_dato="simulacion_teorica" if cat not in ("parametros",) else "parametros_con_fuente",
                            pq="no_declarada_nunca"))
    return out


# ------------------------------------------------------------------------------------------------ 2. verificación
def valores_modelo():
    """Cifras del modelo v0.35 que se comparan (calculadas aquí o leídas de las salidas)."""
    import yaml
    v = json.load(open(ROOT / "out/mc/volumes.json")); M = yaml.safe_load(open(ROOT / "materials.yaml"))
    masa = {}
    for n, i in v.items():
        d = M.get(i["material"], {}).get("density") or 0
        masa[i["material"]] = masa.get(i["material"], 0) + i["volume_cm3"] * d / 1e6
    srg = {c["nombre"]: c for c in json.load(open(ROOT / "validation/validacion_srg.json"))["comprobaciones"]} \
        if "comprobaciones" in json.load(open(ROOT / "validation/validacion_srg.json")) else \
        {c["nombre"]: c for c in json.load(open(ROOT / "validation/validacion_srg.json"))["checks"]}
    vac = json.load(open(ROOT / "out/vacuum/report.json"))["presion"]
    rf = json.load(open(ROOT / "out/fem_rf/resultado_nominal.json"))
    ext = json.load(open(ROOT / "out/fuente_blanco/ext_p0_final.json"))
    gap = json.load(open(ROOT / "validation/gap_efectivo_camara.json"))
    return dict(masa_t=masa, srg=srg, vac=vac, rf=rf, ext=ext, gap=gap)


def fila(id_, parametro, ref, ref_valor, ref_etiqueta, modelo_valor, unidad, fuente_modelo, tol, tol_tipo="abs",
         referencia_doc="", nota="", comparable=True):
    d = None if (modelo_valor is None or ref_valor is None) else modelo_valor - ref_valor
    dr = None if d is None or not ref_valor else d / ref_valor
    if not comparable or d is None:
        coincide = None
    elif tol_tipo == "max":                         # el valor de referencia es una cota superior (p. ej. potencia disponible)
        coincide = modelo_valor <= ref_valor
    else:
        coincide = abs(d) <= (tol if tol_tipo == "abs" else tol * abs(ref_valor))
    cod = str(fuente_modelo).split()[0].split("/")[0].strip("()").upper() if fuente_modelo else ""
    es_iba = cod in FUENTES_IBA
    supuesto = str(ref_etiqueta).startswith("supuesto_declarado")
    # solo un supuesto declarado puede subir de etiqueta; una referencia ya oficial solo confirma o contradice al modelo
    nueva = ETQ_IBA if (coincide and es_iba and supuesto) else ref_etiqueta
    return dict(id=id_, parametro=parametro, referencia=ref, valor_referencia=ref_valor, etiqueta_referencia=ref_etiqueta,
                valor_modelo=modelo_valor, unidad=unidad, fuente_en_modelo=fuente_modelo,
                fuente_IBA=FUENTES_IBA.get(cod), tolerancia=("modelo ≤ referencia" if tol_tipo == "max" else f"±{tol}{'' if tol_tipo == 'abs' else ' (relativa)'}"),
                delta=None if d is None else round(d, 6), delta_rel=None if dr is None else round(dr, 5),
                veredicto=("no_comparable" if coincide is None else "coincide" if coincide else "discrepancia"),
                etiqueta_propuesta=nueva, re_etiquetar=bool(coincide and es_iba and supuesto), documento_oficial=referencia_doc, nota=nota)


def verificacion():
    raw = P_(); V = valores_modelo()
    s = lambda r: cota(raw, r)["src"]
    v = lambda r: cota(raw, r)["v"]
    F = []
    # --- KIUBE-V04 (supuestos de la cámara de vacío de F66/F67/F74) ---
    F.append(fila("V04-1", "radio interior de la cámara de vacío", "KIUBE-V04 (F66/F67/F74)", 548.0, "supuesto_declarado",
                  v("vacuum_chamber.inner_radius"), "mm", s("vacuum_chamber.inner_radius"), 3.0,
                  nota="T04F2 = medida a escala en TUP04 Fig. 2 (±3 mm)"))
    g = V["gap"]
    F.append(fila("V04-2", "gap efectivo de la cámara (volumen de vacío entre caras de polo / π r²)", "KIUBE-V04", 132.0,
                  "supuesto_declarado", round(g["gap_efectivo_mm"], 1), "mm", "CALC (CAD v0.35: colina 24 mm T04F2, valle 200 mm PAT16)", 5.0,
                  nota=("V04 equivale a 0.5·24 + 0.5·240 mm (valle de 240 mm, anterior a v0.21). El modelo: valle Gv = 200 mm "
                        "(PAT16, reiv. 19) → 0.5·24 + 0.5·200 = 112 mm geométrico; con dees, stems y cuñas del corrector en el "
                        f"volumen, {g['gap_efectivo_mm']:.1f} mm (|z| ≤ 100, r < 548). Por anillos: 60–200 mm {g['gap_ef_60_200_mm']:.0f}, "
                        f"200–400 mm {g['gap_ef_200_400_mm']:.0f}, 400–548 mm {g['gap_ef_400_548_mm']:.0f} mm. Volumen entre polos "
                        f"{g['V_camara_r548_m3']:.4f} m³ (V04: {math.pi * 0.548 ** 2 * 0.132:.4f} m³); vacío total con agujeros de bomba y "
                        f"carcasas RF {g['V_vacio_total_m3']:.3f} m³")))
    F.append(fila("V04-3", "número de bocas de bombeo", "KIUBE-V04", 4, "supuesto_declarado", 4, "-",
                  "TUD03 (4 ODP; KIUBE-OM Tabla 1-1: KIUBE 180 = 4 ODP)", 0, referencia_doc="KIUBE-OM Tabla 1-1"))
    F.append(fila("V04-4", "radio de la boca de bombeo", "KIUBE-V04", 127.0, "supuesto_declarado", v("yoke.pump_hole_radius"), "mm",
                  s("yoke.pump_hole_radius"), 3.0))
    F.append(fila("V04-5", "radio de posición de las bocas", "KIUBE-V04", 375.0, "supuesto_declarado", v("yoke.pump_hole_r_pos"), "mm",
                  s("yoke.pump_hole_r_pos"), 3.0))
    F.append(fila("V04-6", "ubicación de las bombas", "KIUBE-V04", None, "supuesto_declarado (EST: 4 en la tapa superior)", None, "-",
                  "DSR + TUD03 + EST", 0, comparable=False,
                  nota=("DISCREPANCIA CUALITATIVA. V04: 4 bocas de bombeo en la tapa superior. Modelo: los 4 agujeros de T04F2 "
                        "atraviesan la tapa A y el disco C en los valles, pero las bombas activas son 2 difusoras bajo la máquina "
                        "(renders de la ficha) y 2 en las carcasas RF (ubicación EST). Cambia la distribución de presión.")))
    F.append(fila("V04-7", "velocidad efectiva de bombeo (H₂O)", "F66/F67/F74 (con V04)", 7559.4, "simulacion_teorica",
                  round(V["vac"]["S_eff_H2O_l_s"], 1), "l/s", "CALC (TPMC del modelo; bomba 6\" EST 1300 l/s H₂O en la brida)",
                  0.10, "rel", nota="V04: s = 0.25 en 4 bocas Ø254 en la tapa; modelo: 4 bombas con baffle/válvula en sus bridas"))
    F.append(fila("V04-8", "presión de equilibrio (H₂O, sin carga de la fuente)", "F74 (con V04)", 1.33e-7, "simulacion_teorica",
                  float(V["vac"]["P_H2O_10h_mbar"]), "mbar", "CALC (TPMC a 10 h de bombeo)", 0.30, "rel"))
    F.append(fila("V04-9", "presión base medida con RF (referencia IBA)", "TUD03T diap. 8 (IBA, medido)", 2.7e-7, "medicion_publicada_IBA",
                  1.33e-7, "mbar", "F74 (V04) frente a la medida", 0.30, "rel", referencia_doc="TUD03T diap. 8",
                  nota="compara la predicción de F74 con la medida de IBA: F74 subestima ×2; el modelo kiube calibra la desgasificación del acero (×30) con esta medida",
                  comparable=True))
    F.append(fila("V04-10", "carga de H₂ de la fuente", "KIUBE-V04 (F66)", 1.0e-3, "supuesto_declarado", float(V["vac"]["Q_h2_mbar_l_s"]),
                  "mbar·l/s", "CALC (calibrada con la transmisión medida 81 %, TUD03)", 0.5, "rel",
                  nota="34 veces la de V04: es lo que explica la pérdida por gas medida (TUD03: 81 % de transmisión)"))
    F.append(fila("V04-11", "pérdida de haz H⁻ por gas residual", "F65/F66 (con V04)", 1.23e-5, "simulacion_teorica",
                  float(raw["vacuum"]["gas_loss_fraction"]["v"]), "fracción", "TUD03 (transmisión fuente→stripper 81 %)", 0.5, "rel",
                  referencia_doc="TUD03 / TUD03T diap. 8",
                  nota="F65/F66: 0.00123 %; IBA mide 19 % de pérdida total y el modelo atribuye 15 % al gas: cuatro órdenes de magnitud"))
    # --- inventario F75 ---
    F.append(fila("F75-1", "masa del acero magnético (yugo + polos + sectores + insertos + tapón)", "F75 (CAD kiube v0.8)", 14.6,
                  "simulacion_teorica", round(V["masa_t"]["steel_1010"], 3), "t", "CALC (CAD v0.35 × 7.87 g/cm³)", 0.02, "rel",
                  nota="ambos salen del CAD kiube; la diferencia es la geometría v0.8 → v0.35 (tapa 181, valle 100, shim de isocronía)"))
    F.append(fila("F75-2", "masa de las bobinas (homogeneizadas)", "F75 (CAD kiube v0.8)", 2.12, "simulacion_teorica",
                  round(V["masa_t"]["copper_coil_homog"], 3), "t", "CALC (Cu 2.00 t de IPAC26 + epoxi/agua 0.12 t EST)", 0.02, "rel",
                  nota="el Cu (2.00 t) es de IPAC26; los 0.12 t de epoxi-vidrio y agua son EST"))
    # --- SRG (ya en los canónicos, rev. F; el modelo usa rev. G) ---
    sg = V["srg"]
    for k, (nombre, ref, tol, uni) in enumerate([("masa del ciclotrón", 18000.0, 0.10, "kg"),
                                                   ("altura total (cima del cáncamo)", 1850.0, 5.0, "mm")]):
        c = sg[nombre]
        F.append(fila(f"SRG-{k + 1}", nombre, "KIUBE-SRG Tabla 1 (canónico RPOS: rev. F 2019; modelo: rev. G 2023)", ref,
                      "fuente_oficial", c["modelo"], uni, "SRG", tol, "rel" if uni == "kg" else "abs",
                      referencia_doc="KIUBE-SRG Tabla 1", nota=c.get("nota", "")))
    for eje in ("X", "Y"):
        c = sg[f"envolvente {eje} (sin servicios ni blancos)"]
        F.append(fila(f"SRG-env{eje}", f"envolvente {eje}", "KIUBE-SRG Tabla 1", 2000.0, "fuente_oficial", c["modelo"], "mm",
                      "SRG", 75.0, referencia_doc="KIUBE-SRG Tabla 1 (2000) y SRG-0080 (máx. 2050)",
                      nota="la ficha DS da 1900 mm; el modelo queda entre ambas fuentes IBA (validar_srg usa 1975 ± 75)"))
    # --- KIUBE-OM (Production Guide, transcrito en docs/fuentes_oficiales_kiube.md) ---
    rf = V["rf"]; ac = rf.get("acoplador", {})
    F.append(fila("OM-1", "bombas ODP (KIUBE 180)", "KIUBE-OM Tabla 1-1", 4, "fuente_oficial", 4, "-", "TUD03", 0,
                  referencia_doc="KIUBE-OM Tabla 1-1"))
    F.append(fila("OM-2", "puertos de salida", "KIUBE-OM cap. 2 / SRG 2.1.x", 8, "fuente_oficial", len(raw["extraction"]["port_offsets"]["v"]) * 4,
                  "-", "TUP04", 0))
    F.append(fila("OM-3", "potencia del amplificador RF frente a la potencia que pide el modelo (paredes con Q real + haz)",
                  "KIUBE-OM Tabla 1-1 (KIUBE 180: 16 kW)", 16.0, "fuente_oficial", round(ac.get("P_total_kW", float("nan")), 2), "kW",
                  "CALC (fem_rf: Q de cobre ideal × factor_Q_real 0.65 EST + haz 180 µA × 18 MeV)", 0.0, "max",
                  referencia_doc="KIUBE-OM Tabla 1-1",
                  nota=(f"con Q ideal: {rf['modos'][0]['P_kW_40kV']:.1f} kW + {ac.get('P_haz_kW', 0):.2f} kW de haz = "
                        f"{rf['modos'][0]['P_kW_40kV'] + ac.get('P_haz_kW', 0):.1f} kW < 16 kW. El amplificador de 16 kW exige un factor "
                        f"de Q real ≥ {rf['modos'][0]['P_kW_40kV'] / (16 - ac.get('P_haz_kW', 0)):.2f} (EST actual 0.65) o un voltaje de dee "
                        "menor que 40 kV: el EST factor_Q_real queda contradicho por la fuente oficial")))
    F.append(fila("OM-4", "energía extraída (en la lámina)", "KIUBE-OM / TUP03: 18 MeV de extracción; 17.5 MeV en el blanco", 18.0,
                  "fuente_oficial", round(V["ext"]["T_MeV"][0], 3), "MeV", "CALC (fem/fuente_a_blanco.py)", 0.25, "abs",
                  nota=(f"dispersión ±{V['ext']['T_MeV'][1]:.2f} MeV. La comparación con 17.5 MeV en el blanco necesita la pérdida en "
                        "las ventanas (12.5 µm + 35 µm de Havar, He): pendiente, no calculada en el modelo")))
    F.append(fila("OM-5", "eficiencia de la lámina stripper", "KIUBE-OM cap. 2 (≈100 %)", 1.0, "fuente_oficial", 1.0, "fracción",
                  "MODELO (cambio de carga sin pérdidas)", 0.01, referencia_doc="KIUBE-OM cap. 2",
                  nota="el modelo la supone, no la calcula: coincidencia por construcción"))
    return F


# ------------------------------------------------------------------------------------------------ 3. impacto en cadena
def impacto(F):
    by = {f["id"]: f for f in F}
    cam = [by[k] for k in ("V04-2", "V04-6", "V04-7", "V04-10") if by[k]["veredicto"] in ("discrepancia", "no_comparable")]
    afect = []
    if cam:
        motivo = "; ".join(f"{f['id']} {f['parametro']}: {f['valor_referencia']} → {f['valor_modelo']} {f['unidad']}"
                           for f in cam)
        afect += [dict(fase="F66", nombre="mapa de presión 3D TPMC", prioridad="alta", motivo=motivo, accion="re-correr con la geometría y el bombeo de kiube v0.35"),
                  dict(fase="F74", nombre="validación analítica Chiggiato", prioridad="alta", motivo=motivo, accion="re-calcular con S_eff, gap y Q de v0.35"),
                  dict(fase="F67", nombre="validación cruzada Molflow+", prioridad="alta", motivo=motivo, accion="regenerar el XML (tools/kiube_camara_molflow.py) o usar out/vacuum/molflow de v0.35"),
                  dict(fase="F65", nombre="pérdida de haz por stripping en gas", prioridad="alta",
                       motivo=f"V04-11: {by['V04-11']['valor_referencia']} frente a {by['V04-11']['valor_modelo']} (TUD03 medido)",
                       accion="usar la carga de H₂ de la fuente y el perfil de presión de v0.35; contrastar con el 81 % medido")]
    afect += [dict(fase="F70", nombre="MC máquina voxelizada + búnker", prioridad="alta",
                   motivo=("voxel de la geometría kiube v0.8. En v0.35: blancos en la cara del yugo (r ≈ 969 mm, antes ≈ 590; VE24), "
                           "puertas de HDPE de 18 cm, acero +2.9 %, valle 100 mm, caja del dee 74.5 mm"),
                   accion="re-voxelizar desde out/mc/kiube_openmc.h5m o volumes.json de v0.35"),
              dict(fase="F71/F72/F73/F76/F77/F78", nombre="derivadas del MC F70 (blindaje, activación, H*(10), mallas, exterior)",
                   prioridad="media", motivo="heredan la geometría de F70", accion="re-correr después de F70"),
              dict(fase="F75", nombre="inventario FEM", prioridad="media",
                   motivo=f"F75-1: {by['F75-1']['delta_rel'] * 100:+.1f} % de acero; 124 volúmenes en v0.35 (137 en v0.8)",
                   accion="regenerar desde out/mc/volumes.json v0.35"),
              dict(fase="F68/F69", nombre="órbitas centrales", prioridad="media",
                   motivo="mapa magnético FE v0.34 (isócrono) y región central con menisco; ΔE/vuelta 100.7 keV (F65: 98.5)",
                   accion="usar out/fem_mag/ngsolve y out/fuente_blanco de v0.35")]
    return afect


# ------------------------------------------------------------------------------------------------ 4. telemetría simulada
def telemetria(dir_):
    raw = P_()
    dir_.mkdir(parents=True, exist_ok=True)
    tp = json.load(open(ROOT / "out/vacuum/tpmc_pumpdown.json"))["pumpdown"]
    with open(dir_ / "bombeo_simulado.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t_h", "P_total_mbar", "P_metales_mbar", "P_juntas_mbar", "origen"])
        for i, t in enumerate(tp["t_h"]):
            w.writerow([t, tp["P_total"][i], tp["P_metales"][i], tp["P_juntas"][i], "simulado_declarado"])
    rf = json.load(open(ROOT / "out/fem_rf/resultado_nominal.json")); fb = json.load(open(ROOT / "validation/fuente_blanco.json"))
    ext = json.load(open(ROOT / "out/fuente_blanco/ext_p0_final.json")); vac = json.load(open(ROOT / "out/vacuum/report.json"))["presion"]
    k = raw["coil"]["escala_iso"]["v"]; NI = raw["coil"]["NI_total"]["v"] * (1 + (k - 1) / 0.24) / 1e3
    Pb = raw["coil"]["power_kW"]["v"] * (NI / (raw["coil"]["NI_total"]["v"] / 1e3)) ** 2
    frac = {r["etapa"]: r["pct"] / 100 for r in fb["fraccion_DC"]}
    f_col = [v for kk, v in frac.items() if kk.startswith("pasa el colimador")][0]
    f_lam = [v for kk, v in frac.items() if kk.startswith("llega a la lámina (con gas")][0]
    estado = dict(NI_kAt=round(NI, 2), P_bobinas_kW=round(Pb, 2), f_RF_MHz=round(rf["modos"][0]["f_MHz"], 4), V_dee_kV=40.0,
                  P_RF_paredes_Qideal_kW=round(rf["modos"][0]["P_kW_40kV"], 2),
                  P_RF_paredes_Qreal_kW=round(rf.get("acoplador", {}).get("P_paredes_kW", float("nan")), 2),
                  P_base_H2O_mbar=vac["P_H2O_10h_mbar"], P_operacion_H2_mbar=vac["P_H2_operacion_mbar"],
                  E_lamina_MeV=ext["T_MeV"][0], sigma_E_MeV=ext["T_MeV"][1], sigma_x_blanco_mm=ext["sigma_x_mm"],
                  sigma_z_blanco_mm=ext["sigma_z_mm"], fraccion_lamina_por_mA_fuente=f_lam, fraccion_colimador_por_mA_fuente=f_col)
    # jornada simulada KIUBE 180, doble haz 2 × 90 µA en blanco (PROP), perfil de operación supuesto (EST)
    I_blanco = 90.0
    I_fuente_req = I_blanco / 1000 / f_col                      # mA por puerto que exige el modelo
    filas = []
    for t in range(0, 181):                                     # minutos: 0–20 RF, 20–30 rampa de haz, 30–150 irradiación, 150–180 fin
        rf_on = t >= 5
        haz = 0.0 if t < 20 else min(1.0, (t - 20) / 10) if t < 150 else max(0.0, 1 - (t - 150) / 5)
        Ib = I_blanco * haz
        filas.append(dict(t_min=t, NI_kAt=estado["NI_kAt"], P_bobinas_kW=estado["P_bobinas_kW"],
                          f_RF_MHz=estado["f_RF_MHz"] if rf_on else 0.0, V_dee_kV=40.0 if rf_on else 0.0,
                          P_RF_kW=(estado["P_RF_paredes_Qreal_kW"] + 2 * Ib * 1e-3 * 18.0) if rf_on else 0.0,
                          P_vacio_mbar=(estado["P_base_H2O_mbar"] + (estado["P_operacion_H2_mbar"] if haz > 0 else 0.0)),
                          I_fuente_mA_por_puerto=round(I_fuente_req * haz, 3),
                          I_lamina_uA_por_puerto=round(I_fuente_req * haz * f_lam * 1e3, 2),
                          I_colimador_uA_por_puerto=round(I_fuente_req * haz * (f_lam - f_col) * 1e3, 2),
                          I_blanco_uA_por_puerto=round(Ib, 2), E_lamina_MeV=round(estado["E_lamina_MeV"], 3) if haz > 0 else "",
                          origen="simulado_declarado"))
    with open(dir_ / "jornada_simulada_kiube180.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)
    meta = dict(
        origen="simulado_declarado", tipo_dato="simulacion_teorica", pq="no_declarada_nunca", modelo=f"kiube {VERSION}",
        advertencia=("Telemetría del gemelo, NO medida de haz real (memoria 3). Puede alimentar el Calibration Layer solo como "
                     "dato simulado declarado (pruebas de la cadena, lado predicho); nunca como eob_medido ni en fit_factor."),
        estado_operacion=estado,
        perfil_operacion="EST: RF a los 5 min, rampa de haz 20–30 min, irradiación hasta 150 min, corte en 5 min (supuesto)",
        escenario="KIUBE 180 doble haz, 2 × 90 µA en blanco (PROP: 180 µA)",
        discrepancia_conocida=(f"Para 90 µA por blanco el modelo pide {I_fuente_req:.1f} mA de fuente por puerto: transmisión al "
                               f"blanco {f_col * 100:.2f} % de la corriente DC (v0.35: pérdidas verticales por νr = 2νz a 14–17 MeV y "
                               "mancha ancha). IBA no publica la corriente de fuente; la cifra no es representativa de la máquina real."),
        canales=dict(t_min="min", NI_kAt="kA·vuelta", P_bobinas_kW="kW", f_RF_MHz="MHz", V_dee_kV="kV", P_RF_kW="kW",
                     P_vacio_mbar="mbar", I_fuente_mA_por_puerto="mA", I_lamina_uA_por_puerto="µA", I_colimador_uA_por_puerto="µA",
                     I_blanco_uA_por_puerto="µA", E_lamina_MeV="MeV"),
        calibration_layer=dict(
            uso_permitido="lado predicho (eob_predicho) y simulacro de la cadena con origen=simulado_declarado",
            uso_prohibido="BatchMeasurement.eob_medido_gbq / fit_factor: requiere medida EOB de un lote real (17.6) — pendiente_entorno",
            reaccion="o18_p_n_f18", equipo="IBA Cyclone KIUBE 180 (gemelo kiube v0.35)",
            parametros_predictor=dict(E_lamina_MeV=estado["E_lamina_MeV"], sigma_E_MeV=estado["sigma_E_MeV"],
                                      I_blanco_uA=I_blanco, t_irradiacion_min=120)))
    json.dump(meta, open(dir_ / "telemetria_simulada_meta.json", "w"), indent=1, ensure_ascii=False, default=float)
    return meta


def accion(r):
    sup = str(r["etiqueta_referencia"]).startswith("supuesto_declarado")
    if r["re_etiquetar"]:
        return "re_etiquetado"
    if r["veredicto"] == "no_comparable":
        return "discrepancia_cualitativa_reportada" if "DISCREPANCIA" in r["nota"] else "sin_cambio_no_comparable"
    if r["veredicto"] == "discrepancia":
        return "discrepancia_reportada" if sup else "discrepancia_modelo_vs_referencia"
    return "sin_cambio_coincidencia_con_cota_no_IBA" if sup else "modelo_concuerda_con_referencia"


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(ROOT / "out/rpos_v035")); a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    meta_tel = telemetria(out / "telemetria_simulada")
    F = verificacion(); I = impacto(F)
    hoy = datetime.datetime.now().isoformat(timespec="seconds")
    json.dump(dict(fecha=hoy, modelo=f"kiube {VERSION}", regla_reetiquetado=(
        f"`{ETQ_IBA}` solo si coincide dentro de la tolerancia declarada Y la cota del modelo tiene fuente IBA; el valor y la "
        "etiqueta anteriores se conservan (append-only); toda diferencia se reporta cuantificada"), filas=F),
        open(out / "verificacion_cruzada.json", "w"), indent=1, ensure_ascii=False, default=float)
    # registro append-only del re-etiquetado (una línea por decisión; no se reescriben líneas existentes)
    reg = out / "registro_reetiquetado.jsonl"
    previas = set()
    if reg.exists():
        previas = {json.loads(l)["clave"] for l in reg.read_text().splitlines() if l.strip()}
    with open(reg, "a") as f:
        for r in F:
            clave = hashlib.md5(json.dumps([r["id"], r["valor_referencia"], r["valor_modelo"], r["etiqueta_propuesta"]],
                                           default=str).encode()).hexdigest()
            if clave in previas:
                continue
            f.write(json.dumps(dict(clave=clave, fecha=hoy, id=r["id"], parametro=r["parametro"],
                                    etiqueta_anterior=r["etiqueta_referencia"], valor_anterior=r["valor_referencia"],
                                    etiqueta_nueva=r["etiqueta_propuesta"] if r["re_etiquetar"] else r["etiqueta_referencia"],
                                    accion=accion(r),
                                    delta=r["delta"], delta_rel=r["delta_rel"], fuente_IBA=r["fuente_IBA"]),
                               ensure_ascii=False, default=float) + "\n")
    json.dump(dict(fecha=hoy, afectadas=I), open(out / "impacto_cadena.json", "w"), indent=1, ensure_ascii=False)
    # manifiesto (al final: incluye los archivos del propio paquete, salvo él mismo)
    arts = artefactos()
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name not in ("manifiesto.json", "MD5SUMS"):
            arts.append(dict(ruta=str(p.relative_to(ROOT)), md5=md5(p), bytes=p.stat().st_size, formato=p.suffix.lstrip(".").upper(),
                             categoria="telemetria_simulada" if "telemetria" in str(p) else "paquete_rpos", unidades="ver meta",
                             generador="herramientas/paquete_rpos.py", modelo="kiube", version_modelo=VERSION,
                             fecha_modificacion=hoy, tipo_dato="simulado_declarado" if "telemetria" in str(p) else "simulacion_teorica",
                             pq="no_declarada_nunca"))
    cats = {}
    for x in arts:
        cats[x["categoria"]] = cats.get(x["categoria"], 0) + 1
    man = dict(fecha=hoy, modelo="kiube", version=VERSION, raiz="(raíz del repositorio kiube)", n_artefactos=len(arts),
               por_categoria=cats, bytes_total=sum(x["bytes"] for x in arts), artefactos=arts)
    json.dump(man, open(out / "manifiesto.json", "w"), indent=1, ensure_ascii=False)
    with open(out / "MD5SUMS", "w") as f:
        for x in arts:
            f.write(f"{x['md5']}  {x['ruta']}\n")
    v = {k: sum(1 for r in F if r["veredicto"] == k) for k in ("coincide", "discrepancia", "no_comparable")}
    print(json.dumps(dict(artefactos=len(arts), por_categoria=cats, verificacion=v,
                          re_etiquetadas=[r["id"] for r in F if r["re_etiquetar"]],
                          discrepancias=[(r["id"], r["delta_rel"]) for r in F if r["veredicto"] == "discrepancia"],
                          afectadas=[(x["fase"], x["prioridad"]) for x in I]), ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
