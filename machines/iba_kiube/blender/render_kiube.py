"""Render fotorrealista Blender Cycles del Cyclone KIUBE (modelo paramétrico público).

Entrada: out/render/scene.json + parts/*.stl (python build.py --lod render).
Estudio blanco (ciclorama + 3 luces de área), materiales PBR por pieza, chaflán gris de las
aristas del yugo por normal de cara. Sin logotipos ni marcas.

    blender -b -P blender/render_kiube.py -- --escena out/render/scene.json \
        --salida out/render/kiube_4k.png --res 3840x2160 --muestras 128 [--vista abierta]

Vistas: exterior (por defecto) | abierta (mitad superior elevada por el sistema de elevación).
REGLA: el render es comunicación visual; la geometría es la del modelo (params.yaml).
"""
from __future__ import annotations

import argparse, json, math, sys, time
from pathlib import Path

import bpy  # type: ignore
from mathutils import Vector  # type: ignore


def _args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    p = argparse.ArgumentParser()
    p.add_argument("--escena", required=True)
    p.add_argument("--salida", required=True)
    p.add_argument("--res", default="3840x2160")
    p.add_argument("--muestras", type=int, default=128)
    p.add_argument("--vista", default="exterior", choices=["exterior", "abierta", "servicios"])
    p.add_argument("--elevacion_mm", type=float, default=650.0)
    p.add_argument("--blend", default=None)
    p.add_argument("--tiles", type=int, default=1, help="N×N regiones (8K con poca RAM); unir con blender/unir_tiles.py")
    p.add_argument("--solape", type=int, default=128, help="solape entre regiones en px")
    p.add_argument("--gpu", default="no", choices=["no", "auto", "OPTIX", "CUDA", "HIP", "ONEAPI", "METAL"],
                   help="auto: OptiX > CUDA > HIP > oneAPI > Metal (RTX: OptiX)")
    return p.parse_args(argv)


def _pbsdf(nombre, base, metal=0.0, rough=0.4, coat=0.0, spec=0.5):
    m = bpy.data.materials.new(nombre)
    m.use_nodes = True
    b = m.node_tree.nodes["Principled BSDF"]
    b.inputs["Base Color"].default_value = (*base, 1.0)
    b.inputs["Metallic"].default_value = metal
    b.inputs["Roughness"].default_value = rough
    for k in ("Coat Weight", "Clearcoat"):
        if k in b.inputs:
            b.inputs[k].default_value = coat
            break
    for k in ("Specular IOR Level", "Specular"):
        if k in b.inputs:
            b.inputs[k].default_value = spec
            break
    return m


_POR_DEFECTO = {   # respaldo si scene.json no trae paleta (escenas anteriores a v0.10)
    "pintura_blanca": dict(base=(0.80, 0.81, 0.83), metal=0.0, rough=0.30, coat=0.4),
    "acero_cepillado": dict(base=(0.62, 0.63, 0.65), metal=1.0, rough=0.33),
    "acero_mecanizado": dict(base=(0.58, 0.59, 0.61), metal=1.0, rough=0.22),
    "inox": dict(base=(0.70, 0.71, 0.73), metal=1.0, rough=0.18),
    "cobre": dict(base=(0.85, 0.50, 0.32), metal=1.0, rough=0.25),
    "cobre_bobina": dict(base=(0.55, 0.28, 0.16), metal=0.6, rough=0.45),
    "plastico_negro": dict(base=(0.03, 0.03, 0.035), metal=0.0, rough=0.45),
    "gris_oscuro": dict(base=(0.16, 0.17, 0.18), metal=0.3, rough=0.4),
    "niobio": dict(base=(0.55, 0.57, 0.62), metal=1.0, rough=0.3),
    "aluminio": dict(base=(0.80, 0.81, 0.83), metal=1.0, rough=0.3),
    "hormigon": dict(base=(0.62, 0.61, 0.58), metal=0.0, rough=0.85),
}


def _materiales(paleta=None):
    """Materiales PBR desde la paleta medida en las fichas técnicas (scene.json → 'paleta', generada desde
    palette.yaml por kiube/render_export.py); sin paleta, los valores por defecto de v0.9."""
    tabla = dict(_POR_DEFECTO)
    if paleta:
        tabla.update(paleta["materiales"])
    out = {}
    for k, v in tabla.items():
        out[k] = _pbsdf(k, tuple(v["base"]), v.get("metal", 0.0), v.get("rough", 0.4), coat=v.get("coat", 0.0))
    if "banda_chaflan" not in out:
        out["banda_chaflan"] = out["acero_cepillado"]
    if "acero_inserto" not in out:
        out["acero_inserto"] = out["acero_mecanizado"]
    return out


def _importar(ruta, nombre):
    bpy.ops.wm.stl_import(filepath=str(ruta))
    o = bpy.context.selected_objects[0]
    o.name = o.data.name = nombre
    o.data.transform(__import__("mathutils").Matrix.Scale(0.001, 4))     # mm → m (unidades físicas de Cycles)
    return o


def _suavizar(o, ang=32.0):
    bpy.context.view_layer.objects.active = o
    o.select_set(True)
    try:
        bpy.ops.object.shade_smooth_by_angle(angle=math.radians(ang))
    except Exception:
        bpy.ops.object.shade_smooth()
    o.select_set(False)


def _chaflan_gris(o, mats, clave="banda_chaflan"):
    """Caras inclinadas (chaflán 45° de las aristas) → banda de chaflán de la paleta (gris azulado en la foto real,
    metálico en el render oficial)."""
    o.data.materials.append(mats[clave])
    for pl in o.data.polygons:
        nz = abs(pl.normal.z)
        if 0.35 < nz < 0.92:
            pl.material_index = 1


def _estudio(az_cam=-49.0):
    """Ciclorama y luces orientados a la cámara (azimut en grados; −49° = encuadre exterior original)."""
    sc = bpy.context.scene
    giro = math.radians(az_cam + 49.0)
    # ciclorama: piso + fondo curvo continuo
    bpy.ops.mesh.primitive_plane_add(size=1)
    cyc = bpy.context.active_object
    cyc.name = "ciclorama"
    me = cyc.data
    import bmesh
    bm = bmesh.new()
    W, D, R, Hh = 26.0, 14.0, 4.0, 12.0
    prof = []
    for i in range(0, 11):                      # piso
        prof.append((-D / 2 + i * (D - R) / 10, -1.17))
    for i in range(1, 25):                      # curva
        a = i / 24 * math.pi / 2
        prof.append((D / 2 - R + R * math.sin(a), -1.17 + R * (1 - math.cos(a))))
    for i in range(1, 6):                       # pared
        prof.append((D / 2, -1.17 + R + i * (Hh - R) / 5))
    vs = [[bm.verts.new((-W / 2 + j * W, y, z)) for (y, z) in prof] for j in (0, 1)]
    for k in range(len(prof) - 1):
        bm.faces.new((vs[0][k], vs[1][k], vs[1][k + 1], vs[0][k + 1]))
    bm.to_mesh(me); bm.free()
    cyc.rotation_euler = (0, 0, math.radians(az_cam + 90.0))   # fondo detrás de la máquina, opuesto a la cámara
    m = _pbsdf("ciclorama", (0.90, 0.92, 0.95), 0.0, 0.22, spec=0.4)
    me.materials.append(m)
    bpy.ops.object.shade_smooth()
    # mundo: luz ambiente blanco-azulada como la ficha
    w = bpy.data.worlds.new("mundo"); sc.world = w; w.use_nodes = True
    bg = w.node_tree.nodes["Background"]
    bg.inputs[0].default_value = (0.82, 0.87, 0.97, 1.0)
    bg.inputs[1].default_value = 0.30
    # luces de área: principal, relleno, contra
    for nombre, loc, size, en in (("principal", (3.6, -3.4, 4.6), 3.5, 900.0),
                                  ("relleno", (-4.8, -2.8, 2.4), 4.5, 320.0),
                                  ("contra", (-1.6, 5.0, 4.2), 3.0, 650.0)):
        L = bpy.data.lights.new(nombre, "AREA"); L.size = size; L.energy = en
        lo = bpy.data.objects.new(nombre, L); sc.collection.objects.link(lo)
        loc = (loc[0] * math.cos(giro) - loc[1] * math.sin(giro), loc[0] * math.sin(giro) + loc[1] * math.cos(giro), loc[2])
        lo.location = loc
        lo.rotation_euler = (Vector((0, 0, 0.1)) - Vector(loc)).to_track_quat("-Z", "Y").to_euler()


def _camara(vista):
    sc = bpy.context.scene
    cam = bpy.data.cameras.new("cam"); co = bpy.data.objects.new("cam", cam)
    sc.collection.objects.link(co); sc.camera = co
    if vista == "abierta":
        co.location, mira, cam.lens = (4.3, -4.9, 3.3), Vector((0, 0, -0.05)), 50
    elif vista == "servicios":                     # encuadre de la foto VE24 p.7: cara de servicio +x y esquina 45°
        co.location, mira, cam.lens = (5.75, 2.55, 1.45), Vector((0, 0, -0.35)), 50
    else:
        co.location, mira, cam.lens = (4.7, -5.4, 1.75), Vector((0, 0, -0.28)), 50
    co.rotation_euler = (mira - co.location).to_track_quat("-Z", "Y").to_euler()
    cam.dof.use_dof = True
    cam.dof.focus_distance = (mira - co.location).length
    cam.dof.aperture_fstop = 11.0
    return co


def _gpu(pedido):
    """Activa la GPU en Cycles. Devuelve 'GPU' o 'CPU' (si no hay dispositivo compatible)."""
    if pedido == "no":
        return "CPU"
    prefs = bpy.context.preferences.addons["cycles"].preferences
    tipos = ["OPTIX", "CUDA", "HIP", "ONEAPI", "METAL"] if pedido == "auto" else [pedido]
    for t in tipos:
        try:
            prefs.compute_device_type = t
        except TypeError:
            continue
        prefs.refresh_devices() if hasattr(prefs, "refresh_devices") else prefs.get_devices()
        gpus = [d for d in prefs.devices if d.type == t]
        if gpus:
            for d in prefs.devices:
                d.use = d.type == t
            print(f"[gpu] {t}: " + ", ".join(d.name for d in gpus), flush=True)
            return "GPU"
    print("[gpu] sin GPU compatible → CPU", flush=True)
    return "CPU"


def main():
    a = _args()
    esc = json.load(open(a.escena)); base = Path(a.escena).parent
    bpy.ops.wm.read_factory_settings(use_empty=True)
    sc = bpy.context.scene
    mats = _materiales(esc.get("paleta"))
    t0 = time.time()
    for p in esc["piezas"]:
        o = _importar(base / p["stl"], p["nombre"])
        o.data.materials.append(mats[p["material"]])
        if p["nombre"].startswith("yoke_AB"):
            _chaflan_gris(o, mats, (esc.get("paleta") or {}).get("chaflan", "banda_chaflan"))
        _suavizar(o)
        if a.vista == "abierta" and p["mitad_superior"]:
            o.location.z += a.elevacion_mm / 1000.0
        if a.vista == "abierta" and p["nombre"].startswith(("port_frame", "ion_source_cover")):
            o.hide_render = True                       # marco y tapa de la fuente atraviesan el plano medio
    if a.vista == "abierta":                           # vástagos extendidos de los gatos de elevación
        for o in [o for o in bpy.data.objects if o.name.startswith("lifting_jack_")]:
            bb = [o.matrix_world @ Vector(c) for c in o.bound_box]
            cx = sum(v.x for v in bb) / 8; cy = sum(v.y for v in bb) / 8; z0 = max(v.z for v in bb)
            z1 = z0 + a.elevacion_mm / 1000.0
            bpy.ops.mesh.primitive_cylinder_add(vertices=48, radius=0.032, depth=z1 - z0, location=(cx, cy, (z0 + z1) / 2))
            rod = bpy.context.active_object; rod.name = f"vastago_{o.name}"
            rod.data.materials.append(mats["inox"]); bpy.ops.object.shade_smooth()
    _estudio(23.9 if a.vista == "servicios" else -49.0)
    _camara(a.vista)
    W, H = (int(v) for v in a.res.lower().split("x"))
    r = sc.render
    r.engine = "CYCLES"; r.resolution_x, r.resolution_y, r.resolution_percentage = W, H, 100
    r.image_settings.file_format = "PNG"; r.image_settings.color_depth = "8"; r.image_settings.compression = 30
    try:
        sc.view_settings.view_transform = "AgX"; sc.view_settings.look = "AgX - Medium High Contrast"
    except Exception:
        sc.view_settings.view_transform = "Filmic"
    sc.view_settings.exposure = -0.35
    cy = sc.cycles
    cy.device = _gpu(a.gpu); cy.samples = a.muestras; cy.use_adaptive_sampling = True; cy.adaptive_threshold = 0.02
    cy.use_denoising = True
    try:
        cy.denoiser = "OPENIMAGEDENOISE"
    except Exception:
        pass
    cy.max_bounces, cy.glossy_bounces, cy.diffuse_bounces = 8, 4, 3
    r.use_persistent_data = True
    r.filepath = str(Path(a.salida).resolve())
    if a.blend:
        bpy.ops.wm.save_as_mainfile(filepath=str(Path(a.blend).resolve()))
    t1 = time.time()
    if a.tiles <= 1:
        bpy.ops.render.render(write_still=True)
    else:
        # regiones con solape: cada una cabe en memoria como un render 4K; se unen después con fundido
        n, ov = a.tiles, a.solape
        r.use_border, r.use_crop_to_border = True, True
        meta = {"W": W, "H": H, "n": n, "solape": ov, "tiles": []}
        stem = Path(a.salida).with_suffix("")
        for i in range(n):
            for j in range(n):
                x0 = max(0, i * W // n - ov); x1 = min(W, (i + 1) * W // n + ov)
                y0 = max(0, j * H // n - ov); y1 = min(H, (j + 1) * H // n + ov)
                r.border_min_x, r.border_max_x = x0 / W, x1 / W
                r.border_min_y, r.border_max_y = y0 / H, y1 / H          # Blender: y desde abajo
                f = f"{stem}_t{i}{j}.png"
                r.filepath = str(Path(f).resolve())
                bpy.ops.render.render(write_still=True)
                meta["tiles"].append({"f": f, "x0": x0, "x1": x1, "y0": y0, "y1": y1})
                print(f"[tile] {i}{j} listo ({time.time() - t1:.0f} s)", flush=True)
        json.dump(meta, open(f"{stem}_tiles.json", "w"), indent=1)
    print(f"[render] {a.vista} {W}x{H} {a.muestras} muestras: carga {t1 - t0:.0f} s, render {time.time() - t1:.0f} s → {a.salida}")


if __name__ == "__main__":
    main()
