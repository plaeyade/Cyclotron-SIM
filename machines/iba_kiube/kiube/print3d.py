"""Exportación para impresión 3D (vía Blender) del Cyclone KIUBE.

Del LOD 'print' (sólidos ≥ 0.5 mm) se forman CUERPOS IMPRIMIBLES (uniones booleanas de
piezas que se tocan), se escalan, se teselan, se validan con trimesh y se exportan:
  - stl/<cuerpo>.stl         mm a escala de impresión, una malla estanca por cuerpo
  - kiube_print.glb          escena con nombres y colores por material (importación directa en Blender)
  - manifest.json            escala, colores, validación, espesor mínimo, orientación sugerida

La geometría es la del modelo paramétrico (fuentes en params.yaml); este módulo no añade
nada físico: solo agrupa, escala y verifica imprimibilidad.
"""
import json, math, pathlib, re
import numpy as np
import cadquery as cq
import trimesh

# color de visualización/impresión multicolor por material (sRGB 0–1)
COLOR = {  # v0.9, por material físico: ahora es la variante «didactico» de kiube/palette.py
         "steel_1010": (0.52, 0.55, 0.60), "copper_coil_homog": (0.80, 0.50, 0.20),
         "copper_ofe": (0.78, 0.33, 0.18), "stainless_304L": (0.78, 0.78, 0.80),
         "aluminum_6061": (0.86, 0.86, 0.90), "niobium": (0.45, 0.47, 0.62), "insert": (0.22, 0.42, 0.85),
         "polymer_black": (0.05, 0.05, 0.06), "hdpe": (0.08, 0.08, 0.09)}

# cuerpo imprimible ← regex de piezas; 'orient' = cara que va a la cama
GROUPS = [
    ("yugo_inferior",   r"^(yoke_AB|F_plug|C_sector_\d|D_pole_\d|dummy_dee_\d_\d)_lower$", "steel_1010", "tapa_abajo"),
    ("yugo_superior",   r"^(yoke_AB|F_plug|C_sector_\d|D_pole_\d|dummy_dee_\d_\d)_upper$", "steel_1010", "tapa_abajo"),
    ("inserto_{k}_{h}", r"^E_insert_(\d)_(upper|lower)$", "insert", "tapa_abajo"),
    ("bobina_{h}",      r"^coil_(upper|lower)$", "copper_coil_homog", "plano_abajo"),
    ("camara_vacio",    r"^vacuum_chamber$", "stainless_304L", "plano_abajo"),
    ("dee_{k}_{h}",     r"^(dee|dee_stem)_(\d)_(upper|lower)$", "copper_ofe", "tapa_abajo"),
    ("puente_dees_{h}", r"^dee_bridge_(upper|lower)$", "copper_ofe", "plano_abajo"),
    ("acoplador_rf",    r"^rf_coupler_(loop|feed)$", "copper_ofe", "plano_abajo"),
    ("sintonizador_rf", r"^rf_tuner_(plunger|actuator)$", "copper_ofe", "eje_vertical"),
    ("blanco_T{k}",     r"^T(\d)_(body|collimator|ve_wedge_flange|tubo_haz)$", "niobium", "eje_vertical"),
    ("puerta_blindada_T{k}", r"^T(\d)_puerta_blindada$", "polymer_black", "plano_abajo"),
    ("stripper_{k}",    r"^stripper_(arm|carousel)_(\d)$", "stainless_304L", "plano_abajo"),
    ("fuente_iones_{k}", r"^ion_source_(?:stem_)?(\d)$", "stainless_304L", "eje_horizontal"),
    ("pata_{k}",        r"^leg_(\d)$", "steel_1010", "plano_abajo"),
    ("gato_elevacion_{k}", r"^(lifting_jack|lifting_eye)_(\d)$", "steel_1010", "plano_abajo"),
    ("cancamo_central", r"^lifting_eyebolt_top$", "steel_1010", "plano_abajo"),
    ("gato_aux_{k}",    r"^lifting_aux_(\d_\d)$", "steel_1010", "plano_abajo"),
    ("carcasa_rf_{k}_{h}", r"^rf_housing_(\d)_(upper|lower)$", "stainless_304L", "tapa_abajo"),
    ("brida_ciega_{k}", r"^blank_flange_(\d)_upper$", "stainless_304L", "plano_abajo"),
    ("bomba_difusora_{k}", r"^diffusion_pump_(\d)$", "stainless_304L", "plano_abajo"),
    ("actuador_fuente_{k}", r"^ion_source_(?:actuator|cover)_(\d)$", "stainless_304L", "plano_abajo"),
    ("marco_puerto_{k}", r"^port_frame_(\d)$", "steel_1010", "plano_abajo"),
    ("faldon",          r"^skirt_", "polymer_black", "plano_abajo"),
    ("anillo_tuberias", r"^water_pipe_ring$", "copper_ofe", "plano_abajo"),
    # servicios de la foto real (kiube/auxiliares.py)
    ("cadena_portacables", r"^cable_chain_\d$", "polymer_black", "plano_abajo"),
    ("servicios_agua", r"^(water_header|water_bracket|water_branch|water_valve_handle|hose_blue|hose_black|flow_tx)",
     "brass_cuzn39", "plano_abajo"),
    ("caja_electrica", r"^elec_(box|bracket)", "cabinet_equiv", "plano_abajo"),
    ("manguera_neumatica", r"^pneumatic_coil$", "hose_pu", "plano_abajo"),
]
HALF = {"upper": "sup", "lower": "inf"}

PRINTER = {"fdm": {"min_wall_mm": 0.8, "nozzle_mm": 0.4},
           "resina": {"min_wall_mm": 0.3, "xy_mm": 0.05}}


def group_parts(parts):
    out = {}
    for p in parts:
        for tmpl, rx, mat, orient in GROUPS:
            m = re.match(rx, p["name"])
            if not m:
                continue
            g = m.groups()
            k = next((x for x in g if x and (x.isdigit() or re.fullmatch(r"\d_\d", x))), "")
            h = next((HALF[x] for x in g if x in HALF), "")
            name = tmpl.format(k=k, h=h)
            out.setdefault(name, dict(name=name, mat=mat, orient=orient, parts=[]))["parts"].append(p)
            break
    return out


def fuse(shapes):
    s = shapes[0]
    for t in shapes[1:]:
        s = s.fuse(t)
    return s.clean()


def to_trimesh(shape, scale, tol_print_mm):
    v, f = shape.tessellate(tol_print_mm / scale, 0.15)
    V = np.array([(p.x, p.y, p.z) for p in v]) * scale
    m = trimesh.Trimesh(V, np.array(f), process=True)
    m.merge_vertices()
    trimesh.repair.fix_normals(m)
    return m


def min_thickness(m, n=1500, seed=1):
    """Espesor local por rayos: desde puntos de superficie hacia adentro (−normal)."""
    if m.is_empty:
        return float("nan"), float("nan")
    pts, fi = trimesh.sample.sample_surface(m, n, seed=seed)
    d = -m.face_normals[fi]
    o = pts + d * 1e-4
    loc, ri, _ = m.ray.intersects_location(o, d, multiple_hits=False)
    if not len(ri):
        return float("nan"), float("nan")
    t = np.linalg.norm(loc - o[ri], axis=1)
    t = t[t > 1e-3]
    return float(np.min(t)), float(np.percentile(t, 2))


def export_print(P, E, build, scale=0.1, printer="fdm", outdir="out/print3d", tol_print_mm=0.03, paleta="foto"):
    from . import palette as PAL
    out = pathlib.Path(outdir)
    (out / "stl").mkdir(parents=True, exist_ok=True)
    parts, info = build(P, E, "print")
    groups = group_parts(parts)
    for i, t in enumerate(E["targets"]):
        if f"blanco_T{i+1}" in groups:
            groups[f"blanco_T{i+1}"]["axis_deg"] = t["dir"]
    lim = PRINTER[printer]["min_wall_mm"]
    scene = trimesh.Scene()
    pin = {"radio_mm": 2.0, "alto_mm": 3.0, "holgura_mm": 0.15 if printer == "fdm" else 0.08,
           "radio_posicion_mm": round(scale * (P["yoke"]["cavity_radius"] +
                                               ((P["yoke"]["side"] - P["yoke"]["corner_chamfer"]) / math.sqrt(2)
                                                - P["yoke"]["cavity_radius"]) / 2), 2),
           "azimuts_deg": [45.0, 135.0, 225.0, 315.0],
           "nota": "en la cara de unión del yugo lateral B (z = 0), sobre las diagonales: lejos de puertos y bobina"}
    manifest = dict(escala=scale, pines_registro=pin, escala_txt=f"1:{round(1 / scale)}", unidades="mm (a escala de impresión)",
                    impresora=printer, espesor_min_impresora_mm=lim,
                    origen="centro de la máquina, plano medio z = 0 (antes de orientar en la cama)",
                    nota=("Geometría del modelo paramétrico público (params.yaml, fuentes por cota). "
                          "Piezas < 0.5 mm reales (láminas, ventanas, fluidos) excluidas del LOD print."),
                    cuerpos=[])
    for name, g in sorted(groups.items()):
        shape = fuse([p["shape"] for p in g["parts"]])
        m = to_trimesh(shape, scale, tol_print_mm)
        tmin, t2 = min_thickness(m)
        cad_vol = shape.Volume() * scale ** 3
        fc, nombres, idx = PAL.colores_cuerpo(m, g["parts"], scale, paleta, mat_fisico=g["mat"],
                                              tess=lambda sh: to_trimesh(sh, scale, 4 * tol_print_mm))
        m.visual = trimesh.visual.ColorVisuals(m, face_colors=fc)
        dom = int(np.bincount(idx, weights=m.area_faces).argmax())          # color predominante (por área)
        c = [round(int(x) / 255, 3) for x in fc[int(np.argmax(idx == dom))][:3]]
        m.export(out / "stl" / f"{name}.stl")
        scene.add_geometry(m, node_name=name, geom_name=name)
        bb = m.bounds
        manifest["cuerpos"].append(dict(
            nombre=name, material=g["mat"], color_rgb=c, orientacion=g["orient"],
            piezas=[p["name"] for p in g["parts"]], stl=f"stl/{name}.stl",
            estanca=bool(m.is_watertight), orientacion_consistente=bool(m.is_winding_consistent),
            volumen_mm3=round(float(m.volume), 2), volumen_cad_mm3=round(cad_vol, 2),
            error_volumen_pct=round(100 * (m.volume / cad_vol - 1), 3) if cad_vol else None,
            caja_mm=[round(float(x), 2) for x in (bb[1] - bb[0])],
            espesor_min_mm=round(tmin, 3), espesor_p2_mm=round(t2, 3),
            imprimible=("ok" if t2 >= 1.05 * lim else "limite" if t2 >= 0.95 * lim else "no"),
            eje_deg=g.get("axis_deg"), mitad=("sup" if name.endswith("_sup") or name == "yugo_superior"
                                               else "inf" if name.endswith("_inf") or name == "yugo_inferior" else None),
            triangulos=int(len(m.faces))))
    scene.export(out / "kiube_print.glb")
    json.dump(manifest, open(out / "manifest.json", "w"), indent=1, ensure_ascii=False)
    return manifest
