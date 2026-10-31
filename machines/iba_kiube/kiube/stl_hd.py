"""STL de alta definición del KIUBE para impresión 3D: por pieza y completo como unidad.

De dónde sale el detalle (honesto): los triángulos se generan directamente desde el B-rep del CAD paramétrico
(superficies exactas) con una tolerancia de cuerda y de ángulo; más polígonos = curvas más lisas, NO más detalle
del que tiene el modelo (que es una reconstrucción de literatura pública, con cotas estimadas donde no hay datos).

Pasos por cuerpo imprimible (los 61 de print3d.GROUPS, que cubren las 196 piezas):
  1. unión CAD de sus piezas (OCC) → teselado con flecha ≤ 'flecha_mm' a la escala de salida y ángulo ≤ 'angulo_rad';
  2. validación: estanco, orientación consistente, variedad (manifold3d), auto-intersecciones (MeshLab),
     error de volumen frente al CAD;
  3. remallado isotrópico opcional (MeshLab, preserva aristas vivas) con su desviación de Hausdorff medida;
  4. STL binario por cuerpo + manifest compatible con blender/blender_prep_impresion.py.
Completo: unión booleana real de los 61 cuerpos (manifold3d) → un solo sólido; además la versión multicuerpo
(sin unir), un 3MF con color (paleta de las fichas técnicas, kiube/palette.py) y un GLB con color por cara.

Color (v0.10): cada cara toma el material visual de la pieza CAD de la que procede (paleta medida en las fichas,
palette.yaml). Un cuerpo con piezas de varios colores (yugo + polos, blanco + colimador, actuador + tapa…) se
exporta además partido por color en sólidos estancos (por_color/), que en el 3MF son las PARTES de un mismo objeto:
el slicer multicolor asigna un filamento a cada parte. La banda de los chaflanes del yugo es color de CARA
(no hay sólido que separar): se ve en visores; para imprimirla, pintar en el slicer o a mano.
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
import trimesh

from .print3d import group_parts, fuse, min_thickness, PRINTER
from . import palette as PAL

# flecha máxima [mm a la escala de salida] y ángulo máximo entre normales de facetas vecinas [rad]
CALIDAD = {"media": (0.02, 0.10), "alta": (0.005, 0.05), "ultra": (0.002, 0.02)}
# ángulo mínimo por cuerpo [rad]: en tubos finos (Ø 1–2 mm a 1:10) el criterio de ángulo, no la flecha, fija el número
# de triángulos (la espiral neumática daba 1.6 M con 0.05 rad). Con estos valores la flecha sigue ≤ la pedida.
ANGULO_MIN = {"manguera_neumatica": 0.30, "servicios_agua": 0.15, "cancamo_central": 0.15}


def tessellate(shape, scale, flecha_mm, angulo_rad):
    """Teselado con flecha ABSOLUTA (OCC BRepMesh con isRelative = False). Ojo: Shape.tessellate de CadQuery usa
    tolerancia RELATIVA al tamaño de cada arista; aquí se malla antes en absoluto y CadQuery reutiliza esa malla."""
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.BRepTools import BRepTools
    defl = flecha_mm / scale                                      # en unidades del modelo (mm reales)
    BRepTools.Clean_s(shape.wrapped)
    BRepMesh_IncrementalMesh(shape.wrapped, defl, False, angulo_rad, True)
    v, f = shape.tessellate(defl, angulo_rad)
    V = np.array([(p.x, p.y, p.z) for p in v]) * scale
    m = trimesh.Trimesh(V, np.array(f), process=True)
    m.merge_vertices()
    trimesh.repair.fix_normals(m)
    return m


def tessellate_group(shapes, scale, flecha_mm, angulo_rad):
    """Robusto: tesela cada SÓLIDO por separado y los une con manifold3d (booleanos exactos). Evita las
    auto-intersecciones que deja la unión CAD de caras coincidentes/tangentes al teselar."""
    import manifold3d as mf
    Ms = []
    for sh in shapes:
        for so in (sh.Solids() if hasattr(sh, "Solids") else [sh]):
            m = tessellate(so, scale, flecha_mm, angulo_rad)
            try:
                Ms.append(to_manifold(m))
            except Exception:                                     # noqa: BLE001
                m2 = repair(m)
                Ms.append(to_manifold(m2))
    U = Ms[0] if len(Ms) == 1 else mf.Manifold.batch_boolean(Ms, mf.OpType.Add)
    u = from_manifold(U)
    trimesh.repair.fix_normals(u)
    return u


def repair(m):
    """Reparación con MeshLab: duplicados, no-variedad, agujeros pequeños."""
    import pymeshlab as ml
    ms = ml.MeshSet(); ms.add_mesh(ml.Mesh(m.vertices, m.faces))
    ms.meshing_remove_duplicate_vertices(); ms.meshing_remove_duplicate_faces()
    ms.meshing_repair_non_manifold_edges(); ms.meshing_close_holes(maxholesize=1000)
    r = ms.current_mesh()
    out = trimesh.Trimesh(r.vertex_matrix(), r.face_matrix(), process=True); trimesh.repair.fix_normals(out)
    return out


def fix_self_intersections(m, dilate=1):
    """MeshLab: selecciona caras auto-intersectadas (+1 anillo), las borra y cierra los agujeros. Devuelve la
    malla y la desviación de Hausdorff respecto de la original."""
    import pymeshlab as ml
    ms = ml.MeshSet(); ms.add_mesh(ml.Mesh(m.vertices, m.faces))
    ms.compute_selection_by_self_intersections_per_face()
    for _ in range(dilate):
        ms.apply_selection_dilatation()
    ms.meshing_remove_selected_faces(); ms.meshing_remove_unreferenced_vertices()
    ms.meshing_close_holes(maxholesize=100000, newfaceselected=False, selfintersection=False)
    r = ms.current_mesh()
    out = trimesh.Trimesh(r.vertex_matrix(), r.face_matrix(), process=True); trimesh.repair.fix_normals(out)
    ms.add_mesh(ml.Mesh(m.vertices, m.faces))
    h = ms.get_hausdorff_distance(sampledmesh=0, targetmesh=1, samplenum=20000)
    return out, dict(hausdorff_max_mm=float(h["max"]), hausdorff_media_mm=float(h["mean"]))


def add_print_bridges(U, radio_mm=1.2, n_muestras=4000):
    """Une al bloque principal cada componente suelto con un puntal cilíndrico entre los puntos más cercanos
    (gatos de elevación, anillo de tuberías, carruseles de stripper: en la máquina cuelgan de soportes que el modelo
    no incluye). Son ayudas de impresión, NO geometría de la máquina: se declaran en el manifest."""
    import manifold3d as mf
    comps = sorted(U.decompose(), key=lambda c: -c.volume())
    # v0.27: descarta láminas degeneradas (espesor ~0, de caras coplanares entre mitades sup/inf): no son material
    def _espesor(c):
        b = from_manifold(c).bounds
        return float((b[1] - b[0]).min())
    comps = [comps[0]] + [c for c in comps[1:] if _espesor(c) > 0.05 and c.volume() > 1e-3]
    if len(comps) < 2:
        return comps[0], []
    main = from_manifold(comps[0])
    pq = trimesh.proximity.ProximityQuery(main)
    out, info = [comps[0]], []
    for c in comps[1:]:
        t = from_manifold(c)
        pts, _ = trimesh.sample.sample_surface(t, n_muestras, seed=2)
        cp, dist, _ = pq.on_surface(pts)
        i = int(np.argmin(dist)); a, b = pts[i], cp[i]
        d = b - a; L = float(np.linalg.norm(d))
        u = d / L if L > 1e-9 else np.array([0.0, 0.0, 1.0])
        p0, p1 = a - 0.6 * radio_mm * u, b + 0.6 * radio_mm * u          # penetra en ambos cuerpos
        cyl = trimesh.creation.cylinder(radius=radio_mm, segment=[p0, p1], sections=32)
        out += [c, to_manifold(cyl)]
        bb = t.bounds
        info.append(dict(componente_caja_mm=[round(float(x), 1) for x in bb[1] - bb[0]],
                         componente_centro_mm=[round(float(x), 1) for x in t.centroid], separacion_mm=round(L, 3),
                         radio_puntal_mm=radio_mm))
    return mf.Manifold.batch_boolean(out, mf.OpType.Add), info


def to_manifold(m):
    import manifold3d as mf
    mesh = mf.Mesh(vert_properties=np.asarray(m.vertices, np.float32), tri_verts=np.asarray(m.faces, np.uint32))
    return mf.Manifold(mesh)


def from_manifold(M):
    # sin fusionar vértices: manifold3d ya entrega la malla indexada y duplica a propósito los vértices donde dos
    # sólidos sólo se tocan en una arista (fusionarlos crearía aristas no-variedad)
    out = M.to_mesh()
    return trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=False)


def self_intersections(m):
    """Auto-intersecciones REALES. MeshLab marca candidatas; las coplanares se verifican con una prueba exacta
    (área de la unión de los triángulos del plano = suma de áreas → no se solapan: falso positivo de triángulos
    finos). Las candidatas no coplanares se cuentan como reales (criterio conservador)."""
    import pymeshlab as ml
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    ms = ml.MeshSet(); ms.add_mesh(ml.Mesh(m.vertices, m.faces))
    ms.compute_selection_by_self_intersections_per_face()
    sel = ms.current_mesh().face_selection_array().astype(bool)
    if not sel.any():
        return 0
    n = m.face_normals; c = m.triangles_center
    key = np.round(np.c_[n, (n * c).sum(1)], 4)
    real = 0
    for k in np.unique(key[sel], axis=0):
        plane = (key == k).all(1)
        if not np.isfinite(k).all() or np.linalg.norm(k[:3]) < 0.5:   # caras degeneradas (área ~0): no solapan
            continue
        nn = k[:3] / np.linalg.norm(k[:3])
        u = np.cross(nn, [0.0, 0.0, 1.0])
        u = u / np.linalg.norm(u) if np.linalg.norm(u) > 1e-6 else np.array([1.0, 0.0, 0.0])
        v = np.cross(nn, u)
        polys = [Polygon(np.c_[t @ u, t @ v]) for t in m.vertices[m.faces[plane]]]
        tot = sum(p.area for p in polys)
        if tot - unary_union(polys).area > 1e-9 * max(tot, 1e-12):
            real += int((plane & sel).sum())
    return real


def remesh_isotropic(m, edge_mm, feature_deg=30.0, iters=6):
    """Remallado isotrópico (MeshLab) preservando aristas vivas; devuelve malla y desviación de Hausdorff [mm]."""
    import pymeshlab as ml
    ms = ml.MeshSet(); ms.add_mesh(ml.Mesh(m.vertices, m.faces))
    ms.meshing_isotropic_explicit_remeshing(iterations=iters, targetlen=ml.PureValue(edge_mm),
                                            featuredeg=feature_deg, checksurfdist=True,
                                            maxsurfdist=ml.PureValue(edge_mm / 10.0))
    r = ms.current_mesh()
    out = trimesh.Trimesh(r.vertex_matrix(), r.face_matrix(), process=True)
    trimesh.repair.fix_normals(out)
    ms.add_mesh(ml.Mesh(m.vertices, m.faces))                    # original como referencia
    h = ms.get_hausdorff_distance(sampledmesh=0, targetmesh=1, samplenum=max(20000, len(out.faces)))
    return out, dict(hausdorff_max_mm=float(h["max"]), hausdorff_media_mm=float(h["mean"]))


def check(m):
    d = dict(triangulos=int(len(m.faces)), estanca=bool(m.is_watertight),
             orientacion_consistente=bool(m.is_winding_consistent))
    try:
        to_manifold(m); d["manifold"] = True
    except Exception:                                             # noqa: BLE001
        d["manifold"] = False
    return d


def export_hd(P, E, build, escala=0.1, calidad="alta", remesh_mm=None, printer="fdm", outdir="out/stl_hd",
              completo=True, autointersecciones=True, puente_mm=1.2, paleta="foto", log=print):
    flecha, ang = CALIDAD[calidad]
    out = pathlib.Path(outdir)
    (out / "piezas").mkdir(parents=True, exist_ok=True)
    parts, _ = build(P, E, "print")
    groups = group_parts(parts)
    for i, t in enumerate(E["targets"]):
        if f"blanco_T{i+1}" in groups:
            groups[f"blanco_T{i+1}"]["axis_deg"] = t["dir"]
    lim = PRINTER[printer]["min_wall_mm"]
    pin = {"radio_mm": 2.0, "alto_mm": 3.0, "holgura_mm": 0.15 if printer == "fdm" else 0.08,
           "radio_posicion_mm": round(escala * (P["yoke"]["cavity_radius"] +
                                                ((P["yoke"]["side"] - P["yoke"]["corner_chamfer"]) / math.sqrt(2)
                                                 - P["yoke"]["cavity_radius"]) / 2), 2),
           "azimuts_deg": [45.0, 135.0, 225.0, 315.0],
           "nota": "en la cara de unión del yugo lateral B (z = 0), sobre las diagonales"}
    man = dict(escala=escala, escala_txt=f"1:{round(1 / escala)}" if escala < 1 else "1:1", unidades="mm (a escala de salida)",
               calidad=calidad, flecha_mm=flecha, angulo_rad=ang, remallado_mm=remesh_mm, impresora=printer,
               espesor_min_impresora_mm=lim, pines_registro=pin,
               origen="centro de la máquina, plano medio z = 0",
               nota=("Triángulos generados desde el CAD exacto (B-rep) con la flecha y el ángulo indicados. "
                     "El detalle es el del modelo paramétrico (cotas de literatura pública; las EST siguen siendo estimadas)."),
               paleta=dict(variante=paleta, archivo="palette.yaml",
                           nota="colores medidos en las fichas IBA (VE24, KEY22) o estimados; ver palette.yaml"),
               cuerpos=[])
    meshes = {}; partes_color = {}; etiquetas_cara = {}
    (out / "por_color").mkdir(exist_ok=True)
    gruesa = lambda sh: tessellate(sh, escala, max(4 * flecha, 0.02), 0.2)   # noqa: E731  (solo para etiquetar caras)
    t0 = time.time()
    for name, g in sorted(groups.items()):
        shape = fuse([p["shape"] for p in g["parts"]])
        ang_g = max(ang, ANGULO_MIN.get(name, 0.0))
        m = tessellate(shape, escala, flecha, ang_g)                 # unión CAD exacta → teselado
        metodo = "union_CAD"
        if not m.is_watertight:                                     # respaldo: sólidos por separado + unión manifold
            m = tessellate_group([p["shape"] for p in g["parts"]], escala, flecha, ang_g); metodo = "union_manifold"
        cad_vol = shape.Volume() * escala ** 3
        info = check(m)
        rem = None
        if remesh_mm:
            m2, rem = remesh_isotropic(m, remesh_mm)
            c2 = check(m2)
            if c2["estanca"] and c2["manifold"]:
                m, info = m2, c2
            else:
                rem["descartado"] = "el remallado rompió la estanqueidad: se conserva el teselado CAD"
        si = self_intersections(m) if autointersecciones else None
        rep = None
        tmin, t2 = min_thickness(m)
        fc, nombres, idx = PAL.colores_cuerpo(m, g["parts"], escala, paleta, mat_fisico=g["mat"], tess=gruesa)
        m.visual = trimesh.visual.ColorVisuals(m, face_colors=fc)
        m.export(out / "piezas" / f"{name}.stl")
        # partido por color en sólidos estancos (para slicers multicolor)
        vis_p = {}
        for pz in g["parts"]:
            vis_p.setdefault(PAL.material_visual(pz["name"], paleta) if paleta != "didactico" else "didactico", []).append(pz)
        sub = []
        if len(vis_p) == 1:
            sub.append(dict(nombre=name, mesh=m, idx=idx, materiales=nombres))
        else:
            for vm, pzs in sorted(vis_p.items()):
                sm = tessellate(fuse([q["shape"] for q in pzs]), escala, flecha, ang_g)
                sfc, snom, sidx = PAL.colores_cuerpo(sm, pzs, escala, paleta, mat_fisico=g["mat"], tess=gruesa)
                sm.visual = trimesh.visual.ColorVisuals(sm, face_colors=sfc)
                f = out / "por_color" / f"{name}__{vm}.stl"
                sm.export(f)
                sub.append(dict(nombre=f"{name}__{vm}", mesh=sm, idx=sidx, materiales=snom, stl=f"por_color/{f.name}",
                                estanca=bool(sm.is_watertight)))
        partes_color[name] = sub
        etiquetas_cara[name] = (m, idx, nombres)
        area = np.bincount(idx, weights=m.area_faces, minlength=len(nombres))
        colores = [dict(**PAL.info(n, paleta), fraccion_area=round(float(a / area.sum()), 4))
                   for n, a in zip(nombres, area)] if paleta != "didactico" else \
                  [dict(material_visual=nombres[0], srgb=[int(x) for x in fc[0][:3]])]
        meshes[name] = m
        bb = m.bounds
        man["cuerpos"].append(dict(
            nombre=name, material=g["mat"], colores=colores, orientacion=g["orient"],
            piezas=[dict(nombre=p["name"], material_visual=PAL.material_visual(p["name"], paleta)) for p in g["parts"]],
            partes_por_color=[dict(nombre=s["nombre"], stl=s.get("stl", f"piezas/{name}.stl"),
                                   estanca=s.get("estanca", True)) for s in sub],
            stl=f"piezas/{name}.stl", **info, caras_autointersectadas=si,
            volumen_mm3=round(float(m.volume), 3), volumen_cad_mm3=round(cad_vol, 3),
            error_volumen_pct=round(100 * (m.volume / cad_vol - 1), 4) if cad_vol else None,
            caja_mm=[round(float(x), 2) for x in (bb[1] - bb[0])], espesor_min_mm=round(tmin, 3), espesor_p2_mm=round(t2, 3),
            imprimible=("ok" if t2 >= 1.05 * lim else "limite" if t2 >= 0.95 * lim else "no"),
            eje_deg=g.get("axis_deg"), remallado=rem, metodo=metodo, reparacion=rep,
            mitad=("sup" if name.endswith("_sup") or name == "yugo_superior"
                   else "inf" if name.endswith("_inf") or name == "yugo_inferior" else None)))
        log(f"[stl_hd] {name}: {info['triangulos']:,} triángulos, estanca {info['estanca']}, manifold {info['manifold']}, "
            f"autointersecciones {si}, Δvol {man['cuerpos'][-1]['error_volumen_pct']} % ({time.time() - t0:.0f} s)")
    if completo:
        (out / "completo").mkdir(exist_ok=True)
        # multicuerpo: todos los cuerpos en un STL, sin unir (cada uno sigue cerrado)
        multi = trimesh.util.concatenate(list(meshes.values()))
        multi.export(out / "completo" / "kiube_completo_multicuerpo.stl")
        # unión booleana real → un solo sólido para imprimir de una pieza
        import manifold3d as mf
        Ms = [to_manifold(m) for m in meshes.values()]
        U = mf.Manifold.batch_boolean(Ms, mf.OpType.Add)
        U, puentes = add_print_bridges(U, radio_mm=puente_mm) if puente_mm else (U, [])
        u = from_manifold(U)
        trimesh.repair.fix_normals(u)
        u.export(out / "completo" / "kiube_completo.stl")
        cu = check(u)
        man["completo"] = dict(stl="completo/kiube_completo.stl", **cu, volumen_mm3=round(float(u.volume), 2),
                               puentes_impresion=puentes,
                               cuerpos_unidos=len(Ms), componentes=int(len(u.split(only_watertight=False))),
                               caja_mm=[round(float(x), 1) for x in (u.bounds[1] - u.bounds[0])],
                               multicuerpo=dict(stl="completo/kiube_completo_multicuerpo.stl", triangulos=int(len(multi.faces))))
        # 3MF con color: un objeto por cuerpo; cuerpos multicolor → partes por color (núcleo 3MF, basematerials)
        objs = []
        for n, sub in partes_color.items():
            partes = []
            for s in sub:
                mats = [(k, PAL.info(k, paleta)["hex"]) if paleta != "didactico" else (k, "#%02X%02X%02X" % tuple(
                    int(x) for x in s["mesh"].visual.face_colors[0][:3])) for k in s["materiales"]]
                partes.append(dict(nombre=s["nombre"], V=s["mesh"].vertices, F=s["mesh"].faces, idx=s["idx"], materiales=mats))
            objs.append(dict(nombre=n, partes=partes))
        r3 = PAL.escribir_3mf(out / "completo" / "kiube_completo.3mf", objs, titulo=f"KIUBE {man['escala_txt']} ({paleta})")
        man["completo"]["3mf"] = dict(archivo="completo/kiube_completo.3mf", **r3,
                                      nota="colores = basematerials del núcleo 3MF por triángulo; cuerpos multicolor "
                                           "como objeto con partes (un filamento por parte en el slicer)")
        # GLB con color por cara (visores, web, Blender)
        PAL.escena_glb(etiquetas_cara, paleta).export(out / "completo" / "kiube_completo_color.glb")
        man["completo"]["glb"] = "completo/kiube_completo_color.glb"
        log(f"[stl_hd] completo: {cu['triangulos']:,} triángulos, estanca {cu['estanca']}, manifold {cu['manifold']}, "
            f"{man['completo']['componentes']} componentes")
    tot = sum(c["triangulos"] for c in man["cuerpos"])
    man["resumen"] = dict(cuerpos=len(man["cuerpos"]), triangulos_total=tot,
                          estancos=sum(c["estanca"] for c in man["cuerpos"]),
                          manifold=sum(c["manifold"] for c in man["cuerpos"]),
                          con_autointersecciones=sum(1 for c in man["cuerpos"] if c["caras_autointersectadas"]),
                          tiempo_s=round(time.time() - t0, 1))
    import csv
    with open(out / "lista_colores.csv", "w", newline="") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["cuerpo", "material_visual", "hex", "R", "G", "B", "RAL_aprox", "dE76", "filamento", "fraccion_area",
                       "confianza", "stl_por_color"])
        for c in man["cuerpos"]:
            stl_c = {s["nombre"].split("__")[-1]: s["stl"] for s in c["partes_por_color"]}
            for k in c["colores"]:
                wcsv.writerow([c["nombre"], k["material_visual"], k.get("hex", ""), *k["srgb"], k.get("ral_aprox", ""),
                               k.get("ral_dE76", ""), k.get("filamento", ""), k.get("fraccion_area", ""), k.get("confianza", ""),
                               stl_c.get(k["material_visual"], c["stl"])])
    json.dump(man, open(out / "manifest.json", "w"), indent=1, ensure_ascii=False)
    return man
