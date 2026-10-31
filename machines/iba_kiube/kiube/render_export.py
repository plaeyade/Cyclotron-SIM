"""Exporta el LOD 'render' para Blender: STL por pieza (mm) + scene.json con material visual y
bandera de 'mitad superior' (para la vista abierta con el sistema de elevación)."""
import json, pathlib, re
import cadquery as cq

from .palette import VIS, CHAFLAN, material_visual, tabla_blender   # paleta medida en las fichas (palette.yaml)
UPPER = re.compile(r"(_upper$|^cable_chain_|^rf_housing_\d_upper|^blank_flange|^lifting_eye"
                   r"|^vacuum_chamber$|^stripper_|^T\d_)")   # la cámara y los blancos suben con la tapa


def export_render(P, E, build, outdir, tol=0.4, variante="foto"):
    out = pathlib.Path(outdir); (out / "parts").mkdir(parents=True, exist_ok=True)
    parts, info = build(P, E, "render")
    scene = {"unidades": "mm", "piso_z_mm": info["floor_z"], "piezas": [],
             "paleta": {"variante": variante if variante != "didactico" else "foto", "fuente": "palette.yaml",
                        "chaflan": CHAFLAN["render" if variante == "render" else "foto"],
                        "materiales": tabla_blender(variante if variante != "didactico" else "foto")}}
    for p in parts:
        vis = material_visual(p["name"], variante if variante != "didactico" else "foto")
        f = out / "parts" / f"{p['name']}.stl"
        cq.exporters.export(p["shape"], str(f), tolerance=tol, angularTolerance=0.08)
        scene["piezas"].append({"nombre": p["name"], "stl": f"parts/{p['name']}.stl", "material": vis,
                                "mitad_superior": bool(UPPER.search(p["name"]))})
    json.dump(scene, open(out / "scene.json", "w"), indent=1, ensure_ascii=False)
    return scene
