"""Preparación para impresión 3D del Cyclone KIUBE en Blender (≥ 4.2; probado en 5.0.1).

Entrada: manifest.json + stl/ que emite `python build.py --lod print3d` (kiube/print3d.py).
Hace, en este orden:
  1. Escena en mm (unidades métricas, scale_length = 0.001) e importación de cada cuerpo.
  2. Colección «Ensamble» (máquina armada, colores por material) para revisar y renderizar.
  3. Verificación manifold con bmesh (aristas no-manifold, vértices sueltos, caras degeneradas).
  4. Pines de registro entre las dos mitades del yugo (macho en la inferior, hembra con
     holgura en la superior), con booleano Manifold/Exact.
  5. Orientación para imprimir (según 'orientacion' del manifest), apoyo en z = 0, división
     en cuadrantes si un cuerpo excede la cama, y acomodo en placas.
  6. Exporta STL por pieza y por placa, guarda el .blend y un informe JSON.

REGLA (igual que los renders de Radiopharma OS): el script NO cambia la geometría del modelo
salvo lo propio de impresión (pines, cortes por cama) y lo declara en el informe. Si falta una
clave del manifest, aborta.

Uso (headless):
    blender -b -P blender/blender_prep_impresion.py -- \
        --manifest out/print3d/manifest.json --salida out/print3d/blender \
        --cama 220x220x250 [--sin-pines] [--blend kiube_impresion.blend]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import bpy  # type: ignore
import bmesh  # type: ignore
from mathutils import Matrix, Vector  # type: ignore


def _args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--salida", required=True)
    p.add_argument("--cama", default="220x220x250", help="volumen de impresión X×Y×Z en mm")
    p.add_argument("--margen", type=float, default=4.0, help="separación entre piezas en la placa (mm)")
    p.add_argument("--sin-pines", action="store_true")
    p.add_argument("--blend", default="kiube_impresion.blend")
    return p.parse_args(argv)


def _req(d, k):
    if k not in d:
        raise SystemExit(f"[abortado] falta la clave '{k}' en el manifest — no se inventa nada")
    return d[k]


# ------------------------------------------------------------------ escena
def _escena_mm():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    us = bpy.context.scene.unit_settings
    us.system, us.scale_length, us.length_unit = "METRIC", 0.001, "MILLIMETERS"


def _coleccion(nombre):
    c = bpy.data.collections.new(nombre)
    bpy.context.scene.collection.children.link(c)
    return c


def _mover_a(obj, col):
    for c in obj.users_collection:
        c.objects.unlink(obj)
    col.objects.link(obj)


def _material(nombre, rgb):
    m = bpy.data.materials.get(nombre) or bpy.data.materials.new(nombre)
    m.use_nodes = True
    b = m.node_tree.nodes.get("Principled BSDF")
    if b:
        b.inputs["Base Color"].default_value = (*rgb, 1.0)
        b.inputs["Metallic"].default_value = 0.6
        b.inputs["Roughness"].default_value = 0.35
    m.diffuse_color = (*rgb, 1.0)
    return m


def _importar_stl(ruta, nombre):
    bpy.ops.wm.stl_import(filepath=str(ruta))
    obj = bpy.context.selected_objects[0]
    obj.name = obj.data.name = nombre
    return obj


# ------------------------------------------------------------------ verificación
def _limpiar(obj, dist=1e-4):
    """Tras un booleano: fusiona dobles, disuelve degenerados y triangula (la malla que se exporta)."""
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=dist)
    for _ in range(4):                                   # hasta que no queden triángulos de área nula
        bmesh.ops.dissolve_degenerate(bm, edges=bm.edges, dist=dist)
        bmesh.ops.triangulate(bm, faces=bm.faces, quad_method="BEAUTY", ngon_method="BEAUTY")
        if not any(f.calc_area() < 1e-8 for f in bm.faces):
            break
        dist *= 5
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def _chequeo(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.triangulate(bm, faces=bm.faces)            # se evalúa la malla tal como se exportará
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    nm = sum(1 for e in bm.edges if not e.is_manifold)
    sueltos = sum(1 for v in bm.verts if not v.link_edges)
    degen = sum(1 for f in bm.faces if f.calc_area() < 1e-8)
    vol = bm.calc_volume(signed=True)
    bm.free()
    return dict(no_manifold=nm, vertices_sueltos=sueltos, caras_degeneradas=degen, volumen_mm3=round(vol, 2))


# ------------------------------------------------------------------ booleanos
def _solver():
    enum = bpy.types.BooleanModifier.bl_rna.properties["solver"].enum_items.keys()
    return "MANIFOLD" if "MANIFOLD" in enum else "EXACT"


def _booleano(obj, cutter, op):
    md = obj.modifiers.new("kiube_bool", "BOOLEAN")
    md.operation, md.object, md.solver = op, cutter, _solver()
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_apply(modifier=md.name)
    _limpiar(obj)


def _cilindro(nombre, r, z0, z1, x, y, verts=64):
    bpy.ops.mesh.primitive_cylinder_add(vertices=verts, radius=r, depth=z1 - z0, location=(x, y, (z0 + z1) / 2))
    o = bpy.context.active_object
    o.name = nombre
    return o


def _pines(inf, sup, pin, informe):
    r, h, hol = _req(pin, "radio_mm"), _req(pin, "alto_mm"), _req(pin, "holgura_mm")
    rp = _req(pin, "radio_posicion_mm")
    for i, az in enumerate(_req(pin, "azimuts_deg")):
        a = math.radians(az)
        x, y = rp * math.cos(a), rp * math.sin(a)
        macho = _cilindro(f"pin_{i}", r, -0.5, h, x, y)                     # penetra 0.5 mm en la inferior
        _booleano(inf, macho, "UNION")
        bpy.data.objects.remove(macho, do_unlink=True)
        hembra = _cilindro(f"hueco_{i}", r + hol, -0.01, h + 0.4, x, y)     # +0.4 mm de fondo
        _booleano(sup, hembra, "DIFFERENCE")
        bpy.data.objects.remove(hembra, do_unlink=True)
    informe["pines"] = dict(n=len(pin["azimuts_deg"]), diametro_mm=2 * r, alto_mm=h, holgura_radial_mm=hol,
                            radio_posicion_mm=rp, solver=_solver(),
                            nota="pines añadidos SOLO para impresión; no forman parte del modelo físico")


# ------------------------------------------------------------------ orientación y cama
def _orientar(obj, c):
    o = _req(c, "orientacion")
    R = Matrix.Identity(4)
    if o == "tapa_abajo" and c.get("mitad") == "sup":
        R = Matrix.Rotation(math.pi, 4, "X")                               # la tapa A queda sobre la cama
    elif o == "eje_vertical" and c.get("eje_deg") is not None:
        R = Matrix.Rotation(math.radians(90), 4, "Y") @ Matrix.Rotation(-math.radians(c["eje_deg"]), 4, "Z")
    elif o == "eje_horizontal":
        R = Matrix.Rotation(math.radians(90), 4, "X")
    obj.data.transform(R)
    obj.data.update()
    bb = [Vector(v.co) for v in obj.data.vertices]
    mn = Vector((min(v.x for v in bb), min(v.y for v in bb), min(v.z for v in bb)))
    mx = Vector((max(v.x for v in bb), max(v.y for v in bb), max(v.z for v in bb)))
    ctr = (mn + mx) / 2
    obj.data.transform(Matrix.Translation((-ctr.x, -ctr.y, -mn.z)))        # centrado en XY, apoyado en z = 0
    obj.data.update()
    return mx - mn


def _dividir_cuadrantes(obj, cama, col):
    """Si no cabe en la cama, se corta en 4 cuadrantes (booleano INTERSECT con cajas)."""
    piezas = []
    L = 10000.0
    for i, (sx, sy) in enumerate(((1, 1), (-1, 1), (-1, -1), (1, -1))):
        cp = obj.copy(); cp.data = obj.data.copy(); col.objects.link(cp)
        cp.name = f"{obj.name}_Q{i+1}"
        bpy.ops.mesh.primitive_cube_add(size=1, location=(sx * L / 2, sy * L / 2, 0))
        caja = bpy.context.active_object
        caja.scale = (L, L, L)
        bpy.context.view_layer.objects.active = caja
        bpy.ops.object.transform_apply(scale=True)
        _booleano(cp, caja, "INTERSECT")
        bpy.data.objects.remove(caja, do_unlink=True)
        piezas.append(cp)
    bpy.data.objects.remove(obj, do_unlink=True)
    return piezas


def _acomodar(objs, cama, margen):
    """Empaquetado por estantes: devuelve {placa: [objs]} y fija la ubicación de cada objeto."""
    W, D, Hc = cama
    orden = sorted(objs, key=lambda o: -o.dimensions.y * o.dimensions.x)
    placas, x, y, fila, p = {}, margen, margen, 0.0, 0
    for o in orden:
        w, d = o.dimensions.x, o.dimensions.y
        if x + w + margen > W:
            x, y, fila = margen, y + fila + margen, 0.0
        if y + d + margen > D:
            p, x, y, fila = p + 1, margen, margen, 0.0
        o.location = (x + w / 2 - W / 2, y + d / 2 - D / 2, 0.0)
        o["placa"] = p
        placas.setdefault(p, []).append(o)
        x += w + margen
        fila = max(fila, d)
    return placas


def _export_stl(objs, ruta):
    bpy.ops.object.select_all(action="DESELECT")
    for o in objs:
        o.select_set(True)
    bpy.context.view_layer.objects.active = objs[0]
    bpy.ops.wm.stl_export(filepath=str(ruta), export_selected_objects=True, apply_modifiers=True,
                          ascii_format=False, global_scale=1.0)


# ------------------------------------------------------------------ principal
def main():
    a = _args()
    man_path = Path(a.manifest)
    man = json.load(open(man_path))
    base = man_path.parent
    out = Path(a.salida)
    (out / "piezas").mkdir(parents=True, exist_ok=True)
    (out / "placas").mkdir(parents=True, exist_ok=True)
    cama = tuple(float(v) for v in a.cama.lower().split("x"))
    _escena_mm()
    col_ens, col_imp = _coleccion("Ensamble"), _coleccion("Impresion")
    informe = dict(escala=_req(man, "escala_txt"), impresora=_req(man, "impresora"), cama_mm=cama,
                   blender=bpy.app.version_string, cuerpos={}, avisos=[])
    objs = {}
    for c in _req(man, "cuerpos"):
        n = _req(c, "nombre")
        o = _importar_stl(base / _req(c, "stl"), n)
        o.data.materials.append(_material(_req(c, "material"), tuple(_req(c, "color_rgb"))))
        # copia armada para el ensamble
        e = o.copy(); e.data = o.data.copy(); e.name = f"{n}_ensamble"
        col_ens.objects.link(e)
        _mover_a(o, col_imp)
        objs[n] = (o, c)
        informe["cuerpos"][n] = dict(chequeo_importado=_chequeo(o), espesor_p2_mm=c["espesor_p2_mm"],
                                     imprimible=c["imprimible"])
        if c["imprimible"] != "ok":
            informe["avisos"].append(f"{n}: espesor p2 {c['espesor_p2_mm']} mm frente a mínimo "
                                     f"{man['espesor_min_impresora_mm']} mm ({c['imprimible']}) — "
                                     "usar resina, subir escala o reforzar")
    if not a.sin_pines and "yugo_inferior" in objs and "yugo_superior" in objs:
        _pines(objs["yugo_inferior"][0], objs["yugo_superior"][0], _req(man, "pines_registro"), informe)
    finales = []
    for n, (o, c) in objs.items():
        dims = _orientar(o, c)
        if dims.x > cama[0] - 2 * a.margen or dims.y > cama[1] - 2 * a.margen:
            informe["avisos"].append(f"{n}: {dims.x:.0f}×{dims.y:.0f} mm excede la cama → 4 cuadrantes")
            for q in _dividir_cuadrantes(o, cama, col_imp):
                _orientar(q, dict(orientacion="plano_abajo"))
                finales.append(q)
        else:
            finales.append(o)
        if dims.z > cama[2]:
            informe["avisos"].append(f"{n}: altura {dims.z:.0f} mm excede la cama en Z")
    bpy.context.view_layer.update()
    for o in finales:
        informe["cuerpos"].setdefault(o.name, {})["chequeo_final"] = _chequeo(o)
        _export_stl([o], out / "piezas" / f"{o.name}.stl")
    placas = _acomodar(finales, cama, a.margen)
    for p, ol in placas.items():
        _export_stl(ol, out / "placas" / f"placa_{p+1:02d}.stl")
    informe["placas"] = {f"placa_{p+1:02d}": [o.name for o in ol] for p, ol in placas.items()}
    col_ens.hide_viewport = col_ens.hide_render = True
    bpy.ops.wm.save_as_mainfile(filepath=str(out / a.blend))
    malos = {k: v["chequeo_final"] for k, v in informe["cuerpos"].items()
             if "chequeo_final" in v and (v["chequeo_final"]["no_manifold"] or v["chequeo_final"]["caras_degeneradas"])}
    informe["resumen"] = dict(piezas=len(finales), placas=len(placas), piezas_no_manifold=len(malos),
                              estado="listo_para_slicer" if not malos else "revisar_no_manifold")
    json.dump(informe, open(out / "informe_impresion.json", "w"), indent=1, ensure_ascii=False)
    print(json.dumps(informe["resumen"], ensure_ascii=False))
    for w in informe["avisos"]:
        print("  aviso:", w)


if __name__ == "__main__":
    main()
