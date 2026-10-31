"""Paleta de colores del KIUBE medida en las fichas técnicas (palette.yaml, generado por palette_src/muestrear.py)
y asignación pieza → material visual, común a render (Blender), impresión 3D (3MF/GLB/STL por color) y visores.

Variantes:
  foto      (por defecto) máquina real fotografiada en VE24 p.7 / KEY22 p.4;
  render    render oficial de VE24 p.1 (chaflanes y gatos metálicos, marcos gris oscuro, parte superior inox);
  didactico colores por material físico de v0.9 (acero azul-gris, insertos azules…): útil para docencia, NO es el real.
"""
from __future__ import annotations

import functools, pathlib, re
import numpy as np
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
VARIANTES = ("foto", "render", "didactico")

# (regex de nombre de pieza, material visual); el primero que coincide gana
VIS = [
    # servicios de la foto real (kiube/auxiliares.py)
    (r"^dummy_dee", "cobre"), (r"^T\d_puerta_blindada", "plastico_negro"), (r"^T\d_tubo_haz", "inox"), (r"^(rf_coupler_loop|rf_tuner_plunger)", "cobre"), (r"^(rf_coupler_feed|rf_tuner_actuator)", "inox"), (r"^lifting_eyebolt", "acero_cepillado"), (r"^cable_chain", "plastico_negro"), (r"^(water_header|water_branch)", "laton"),
    (r"^water_valve_handle", "valvula_roja"), (r"^(hose_blue|pneumatic_coil)", "manguera_azul"),
    (r"^hose_black", "plastico_negro"), (r"^flow_tx", "transmisor_azul"), (r"^elec_box", "caja_electrica_gris"),
    (r"^(water_bracket|elec_bracket)", "gris_oscuro"),
    (r"^yoke_AB", "pintura_blanca"), (r"^(C_sector|D_pole|F_plug|E_insert)", "acero_mecanizado"),
    (r"^coil", "cobre_bobina"), (r"^(dee|dee_stem)", "cobre"), (r"^vacuum_chamber", "inox"),
    (r"^rf_housing", "gris_azulado_claro"), (r"^(blank_flange|diffusion_pump|ion_source_stem)", "inox"),
    (r"^leg", "acero_cepillado"), (r"^lifting", "gris_azulado"), (r"^skirt", "plastico_negro"),
    (r"^water_pipe", "laton"), (r"^port_frame", "marco_puerto"), (r"^ion_source_cover", "pintura_blanca"),
    (r"^ion_source_actuator", "gris_oscuro"), (r"^T\d_(body)", "niobio"), (r"^T\d_ve_wedge", "inox"),
    (r"^T\d_collimator", "aluminio"), (r"^stripper", "inox"), (r"^ion_source_\d$", "cobre"),
    (r"^self_shield", "hormigon"), (r"^(chimenea|chimney)", "cuw"),
]
# diferencias del render oficial respecto de la foto
VIS_RENDER = {"gris_azulado_claro": "inox", "gris_azulado": "acero_cepillado", "laton": "inox"}
CHAFLAN = {"foto": "banda_chaflan", "render": "acero_cepillado", "didactico": None}   # render oficial: banda metálica
# v0.9 (por material físico) para la variante didáctica
DIDACTICO = {"steel_1010": (0.52, 0.55, 0.60), "copper_coil_homog": (0.80, 0.50, 0.20), "copper_ofe": (0.78, 0.33, 0.18),
             "stainless_304L": (0.78, 0.78, 0.80), "aluminum_6061": (0.86, 0.86, 0.90), "niobium": (0.45, 0.47, 0.62),
             "insert": (0.22, 0.42, 0.85), "polymer_black": (0.05, 0.05, 0.06), "brass_cuzn39": (0.80, 0.68, 0.35),
             "hose_pu": (0.15, 0.45, 0.85), "cabinet_equiv": (0.85, 0.85, 0.85), "hdpe": (0.08, 0.08, 0.09)}


@functools.lru_cache(maxsize=None)
def cargar(path=None):
    return yaml.safe_load(open(path or ROOT / "palette.yaml"))


def material_visual(nombre, variante="foto"):
    vis = next((v for rx, v in VIS if re.match(rx, nombre)), "inox")
    if variante == "render":
        vis = VIS_RENDER.get(vis, vis)
    return vis


def srgb(vis, variante="foto"):
    """Color plano sRGB 0–255 del material visual en la variante pedida (si la variante no tiene muestra propia,
    el color principal de la paleta)."""
    e = cargar()["materiales"][vis]
    v = e.get("variantes", {}).get("render") if variante == "render" else None   # 'foto' = color principal
    return tuple(int(x) for x in (v["srgb"] if v else e["srgb"]))


def info(vis, variante="foto"):
    e = cargar()["materiales"][vis]
    c = srgb(vis, variante)
    return dict(material_visual=vis, hex="#%02X%02X%02X" % c, srgb=list(c), ral_aprox=e["ral_aprox"],
                ral_dE76=e["ral_dE76"], filamento=e["filamento"], confianza=e["confianza"], fuente=e["fuente"])


def pbr(vis, variante="foto"):
    """{base (lineal), metal, rough, coat} para el Principled BSDF de Blender."""
    e = cargar()["materiales"][vis]
    v = (e.get("variantes", {}).get("render") if variante == "render" else None) or {}
    p = e["pbr"]
    base = v.get("base_lineal") or p["base_lineal"]
    return dict(base=[float(x) for x in base], metal=float(p["metalico"]), rough=float(p["rugosidad"]),
                coat=float(p.get("barniz", 0.0)))


def tabla_blender(variante="foto"):
    return {k: pbr(k, variante) for k in cargar()["materiales"]}


# ---------------------------------------------------------------- color por cara de una malla impresa
def es_chaflan(normales):
    nz = np.abs(normales[:, 2])
    return (nz > 0.35) & (nz < 0.92)


def colores_cuerpo(m, piezas, escala, variante="foto", mat_fisico=None, n_muestras=60000, tess=None):
    """Color RGBA por cara de la malla de un cuerpo impreso que agrupa varias piezas CAD.
    Cada cara toma el material de la pieza CAD más cercana (KD-tree sobre muestras de la superficie de cada pieza);
    en el yugo, las caras de los chaflanes de 45° (0.35 < |nz| < 0.92) reciben la banda de chaflán.
    Devuelve (colores Nx4 uint8, etiqueta por cara (lista de materiales visuales), índices)."""
    if variante == "didactico":
        c = [int(255 * x) for x in DIDACTICO.get(mat_fisico, (0.6, 0.6, 0.6))] + [255]
        return np.tile(c, (len(m.faces), 1)).astype(np.uint8), [f"didactico:{mat_fisico}"], np.zeros(len(m.faces), int)
    vis = [material_visual(p["name"], variante) for p in piezas]
    etiquetas = sorted(set(vis))
    if len(etiquetas) == 1:
        idx = np.zeros(len(m.faces), int)
    else:
        from scipy.spatial import cKDTree
        pts, lab = [], []
        for p, v in zip(piezas, vis):
            pm = tess(p["shape"]) if tess else None
            if pm is None or pm.is_empty:
                continue
            import trimesh
            s, _ = trimesh.sample.sample_surface_even(pm, max(2000, n_muestras // len(piezas)), seed=1) \
                if pm.area > 0 else (np.zeros((0, 3)), None)
            s = np.vstack([s, pm.vertices])
            pts.append(s); lab.append(np.full(len(s), etiquetas.index(v)))
        tree = cKDTree(np.vstack(pts)); L = np.concatenate(lab)
        _, j = tree.query(m.triangles_center, k=1)
        idx = L[j]
    names = list(etiquetas)
    ch = CHAFLAN.get(variante)
    if ch and any(re.match(r"^yoke_AB", p["name"]) for p in piezas):
        sel = es_chaflan(m.face_normals) & (np.array(names)[idx] == "pintura_blanca")
        if sel.any():
            if ch not in names:
                names.append(ch)
            idx = idx.copy(); idx[sel] = names.index(ch)
    pal = np.array([list(srgb(n, variante)) + [255] for n in names], np.uint8)
    return pal[idx], names, idx


# ---------------------------------------------------------------- 3MF con color (núcleo 3MF: basematerials)
def escribir_3mf(path, objetos, unidad="millimeter", titulo="KIUBE"):
    """objetos: lista de dict(nombre, partes=[dict(nombre, V, F, idx_mat por cara, materiales=[(nombre, hex)])]).
    Cada objeto con varias partes se escribe como objeto con <components> → los slicers (PrusaSlicer, Bambu Studio,
    OrcaSlicer) lo importan como UNA pieza con varias partes, a las que se asigna un filamento cada una.
    Los colores van como <basematerials> del núcleo 3MF (displaycolor) con material por triángulo (pid/p1):
    los muestran Windows 3D Builder/Visor 3D, Blender (importador 3MF) y la mayoría de visores; los slicers
    multicolor asignan el filamento por parte (el color de la parte se toma del material predominante)."""
    import zipfile, io
    mats, mid = [], {}
    for o in objetos:
        for p in o["partes"]:
            for n, h in p["materiales"]:
                if n not in mid:
                    mid[n] = len(mats); mats.append((n, h))
    out = io.StringIO()
    w = out.write
    w('<?xml version="1.0" encoding="UTF-8"?>\n<model unit="%s" xml:lang="es-ES" '
      'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">\n' % unidad)
    w(f'<metadata name="Title">{titulo}</metadata>\n<metadata name="Application">kiube (CadQuery) — stl_hd</metadata>\n')
    w('<resources>\n<basematerials id="1">\n')
    for n, h in mats:
        w(f'<base name="{n}" displaycolor="{h}FF" />\n')
    w('</basematerials>\n')
    oid = 2; build = []
    for o in objetos:
        ids = []
        for p in o["partes"]:
            loc = [mid[n] for n, _ in p["materiales"]]
            idx = np.asarray(p["idx"], int)
            F = np.asarray(p["F"], int)
            ok = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])   # 3MF prohíbe triángulos degenerados
            F, idx = F[ok], idx[ok]
            dom = int(np.bincount(idx, minlength=len(loc)).argmax())
            nm = p["nombre"].replace('"', "")
            w(f'<object id="{oid}" name="{nm}" type="model" pid="1" pindex="{loc[dom]}">\n<mesh>\n<vertices>\n')
            w("".join(f'<vertex x="{v[0]:.5f}" y="{v[1]:.5f}" z="{v[2]:.5f}" />\n' for v in p["V"]))
            w("</vertices>\n<triangles>\n")
            w("".join(f'<triangle v1="{a}" v2="{b}" v3="{c}" pid="1" p1="{loc[k]}" />\n'
                      for (a, b, c), k in zip(F, idx)))
            w("</triangles>\n</mesh>\n</object>\n")
            ids.append(oid); oid += 1
        if len(ids) == 1:
            build.append(ids[0])
        else:
            w(f'<object id="{oid}" name="{o["nombre"]}" type="model">\n<components>\n')
            w("".join(f'<component objectid="{i}" />\n' for i in ids))
            w("</components>\n</object>\n")
            build.append(oid); oid += 1
    w("</resources>\n<build>\n" + "".join(f'<item objectid="{i}" />\n' for i in build) + "</build>\n</model>\n")
    ct = ('<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml" />'
          '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml" /></Types>')
    rel = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
           '<Relationship Target="/3D/3dmodel.model" Id="rel0" '
           'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" /></Relationships>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", rel)
        z.writestr("3D/3dmodel.model", out.getvalue())
    return dict(materiales=len(mats), objetos=len(build))


def escena_glb(cuerpos, variante="foto"):
    """cuerpos: {nombre: (malla, idx por cara, nombres de material)} → trimesh.Scene con un material PBR glTF por
    color (baseColor lineal, metálico, rugosidad de la paleta) y el eje Z del modelo pasado a +Y (convención glTF).
    Se parte cada cuerpo por material en vez de usar colores de vértice: así no se mezclan colores en las aristas."""
    import trimesh
    from trimesh.visual.material import PBRMaterial
    sc = trimesh.Scene()
    zup_a_yup = trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0])
    mats = {}
    for nombre, (m, idx, nombres) in cuerpos.items():
        for k, vis in enumerate(nombres):
            caras = np.nonzero(idx == k)[0]
            if not len(caras):
                continue
            if vis not in mats:
                if variante == "didactico" or vis.startswith("didactico"):
                    c = list(m.visual.face_colors[caras[0]][:3] / 255.0)
                    mats[vis] = PBRMaterial(name=vis, baseColorFactor=c + [1.0], metallicFactor=0.0, roughnessFactor=0.5)
                else:
                    p = pbr(vis, variante)
                    mats[vis] = PBRMaterial(name=vis, baseColorFactor=p["base"] + [1.0], metallicFactor=p["metal"],
                                            roughnessFactor=p["rough"])
            s = m.submesh([caras], append=True)
            s.visual = trimesh.visual.TextureVisuals(material=mats[vis])
            sc.add_geometry(s, node_name=f"{nombre}__{vis}", geom_name=f"{nombre}__{vis}", transform=zup_a_yup)
    return sc
