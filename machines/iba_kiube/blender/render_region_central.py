"""Render de la región central (electrodos + órbitas calculadas) en Blender Cycles.
  blender -b -P blender/render_region_central.py -- --escena .../escena.json --salida rc.png \
          [--res 3840x2160] [--muestras 128] [--gpu auto] [--vista oblicua|planta]
Las placas superiores de dees/dummy dees y las colinas superiores se ocultan para ver las órbitas."""
import argparse, json, math, sys
from pathlib import Path
import bpy  # type: ignore
from mathutils import Vector, Matrix  # type: ignore


def args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    p = argparse.ArgumentParser()
    p.add_argument("--escena", required=True); p.add_argument("--salida", required=True)
    p.add_argument("--res", default="3840x2160"); p.add_argument("--muestras", type=int, default=128)
    p.add_argument("--gpu", default="no"); p.add_argument("--vista", default="oblicua", choices=["oblicua", "planta"])
    return p.parse_args(argv)


def mat(nombre, col, metal=0.0, rough=0.35, emis=0.0):
    m = bpy.data.materials.new(nombre); m.use_nodes = True
    b = m.node_tree.nodes["Principled BSDF"]
    b.inputs["Base Color"].default_value = (*col, 1); b.inputs["Metallic"].default_value = metal
    b.inputs["Roughness"].default_value = rough
    if emis:
        for k in ("Emission Color", "Emission"):
            if k in b.inputs:
                b.inputs[k].default_value = (*col, 1); break
        b.inputs["Emission Strength"].default_value = emis
    return m


def gpu(pedido):
    if pedido == "no":
        return "CPU"
    pr = bpy.context.preferences.addons["cycles"].preferences
    for t in (["OPTIX", "CUDA", "HIP", "ONEAPI", "METAL"] if pedido == "auto" else [pedido]):
        try:
            pr.compute_device_type = t
        except TypeError:
            continue
        pr.refresh_devices() if hasattr(pr, "refresh_devices") else pr.get_devices()
        g = [d for d in pr.devices if d.type == t]
        if g:
            for d in pr.devices:
                d.use = d.type == t
            print("[gpu]", t, [d.name for d in g]); return "GPU"
    print("[gpu] sin GPU → CPU"); return "CPU"


def main():
    a = args(); E = json.load(open(a.escena)); base = Path(a.escena).parent
    bpy.ops.wm.read_factory_settings(use_empty=True); sc = bpy.context.scene
    M = {"cobre": mat("cobre", (0.86, 0.52, 0.33), 1.0, 0.28), "acero": mat("acero", (0.16, 0.17, 0.19), 0.2, 0.6),
         "acero_oscuro": mat("acero_oscuro", (0.35, 0.37, 0.40), 1.0, 0.45),
         "cuw": mat("cuw", (0.55, 0.42, 0.36), 1.0, 0.35)}   # chimeneas de Cu-W (CAS 2017 §3.3.1)
    for p in E["piezas"]:
        bpy.ops.wm.stl_import(filepath=str(base / p["stl"]))
        o = bpy.context.selected_objects[0]
        o.data.transform(Matrix.Scale(0.001, 4))
        o.data.materials.append(M[p["material"]])
        zc = sum((o.matrix_world @ Vector(c)).z for c in o.bound_box) / 8
        if p["grupo"] in ("dees", "dummy_dees", "colinas_tapon", "pullers") and zc > 0:
            o.hide_render = True                      # mitad superior oculta: se ven las órbitas
        try:
            bpy.ops.object.shade_smooth_by_angle(angle=math.radians(30))
        except Exception:
            pass
    for i, c in enumerate(E["orbitas"]):
        col = tuple(c["color"])
        if c["tipo"] == "marcador":
            x, y, z = c["xyz"][0]
            bpy.ops.mesh.primitive_uv_sphere_add(radius=0.0016, location=(x / 1000, y / 1000, z / 1000 + 0.001))
            bpy.context.active_object.data.materials.append(mat(f"m{i}", col, 0, 0.4, 1.2)); continue
        cu = bpy.data.curves.new(f"orb{i}", "CURVE"); cu.dimensions = "3D"
        sp = cu.splines.new("POLY"); pts = c["xyz"]; sp.points.add(len(pts) - 1)
        for k, (x, y, z) in enumerate(pts):
            sp.points[k].co = (x / 1000, y / 1000, z / 1000, 1)
        if c["tipo"] == "haz" and i % 3:
            continue
        cu.bevel_depth = 0.0008 if c["tipo"] == "diseno" else 0.00018; cu.bevel_resolution = 3
        o = bpy.data.objects.new(f"orb{i}", cu); sc.collection.objects.link(o)
        o.data.materials.append(mat(f"o{i}", col, 0, 0.4, 1.6 if c["tipo"] == "diseno" else 0.5))
    w = bpy.data.worlds.new("w"); sc.world = w; w.use_nodes = True
    w.node_tree.nodes["Background"].inputs[0].default_value = (0.02, 0.022, 0.03, 1)
    w.node_tree.nodes["Background"].inputs[1].default_value = 1.0
    for n_, loc, en in (("k", (0.25, -0.3, 0.45), 14.0), ("f", (-0.3, 0.2, 0.3), 5.0)):
        L = bpy.data.lights.new(n_, "AREA"); L.size = 0.3; L.energy = en
        lo = bpy.data.objects.new(n_, L); sc.collection.objects.link(lo); lo.location = loc
        lo.rotation_euler = (Vector((0, 0, 0)) - Vector(loc)).to_track_quat("-Z", "Y").to_euler()
    cam = bpy.data.cameras.new("c"); co = bpy.data.objects.new("c", cam); sc.collection.objects.link(co); sc.camera = co
    co.location = (0.0, 0.0, 0.42) if a.vista == "planta" else (0.17, -0.26, 0.24)
    co.rotation_euler = (Vector((0, -0.005, 0)) - co.location).to_track_quat("-Z", "Y").to_euler()
    cam.lens = 50
    W, H = (int(v) for v in a.res.split("x"))
    r = sc.render; r.engine = "CYCLES"; r.resolution_x, r.resolution_y = W, H
    sc.cycles.device = gpu(a.gpu); sc.cycles.samples = a.muestras; sc.cycles.use_denoising = True
    try:
        sc.view_settings.view_transform = "Standard"
    except Exception:
        pass
    r.filepath = str(Path(a.salida).resolve()); bpy.ops.render.render(write_still=True)


main()
