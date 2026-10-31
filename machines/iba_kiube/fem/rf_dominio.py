"""Dominio de RF del Cyclone KIUBE: el vacío entre conductores donde vive el modo acelerador.

Mitad superior (z ≥ 0; el plano medio es pared magnética para el modo acelerador):
  cámara (r ≤ R_cámara, 0 ≤ z ≤ zA) ∪ agujeros de los valles a través de C y A (0 ≤ z ≤ H)
  ∪ interior de las carcasas RF sobre los valles de los dees (H ≤ z ≤ H + h_carcasa − 20)
  − conductores: yugo (tapa A), sectores C, polos D, insertos E, tapón F, dees, stems, carcasas, bridas ciegas.
Se omiten fuentes de iones, strippers, colimadores y blancos (perturbaciones pequeñas del modo a 40 MHz).
El stem llega a la tapa de la carcasa: resonador coaxial λ/4 cortocircuitado arriba y cargado por el dee.
"""
from __future__ import annotations

import math, pathlib, sys
import cadquery as cq
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from kiube.geometry import build, zcyl                                           # noqa: E402

CONDUCTORES = ("yoke_AB_upper", "C_sector_", "D_pole_", "E_insert_", "F_plug_upper", "dee_", "dee_stem_",
               "rf_housing_", "blank_flange_", "dummy_dee_", "rf_tuner_plunger")
# el lazo del acoplador (Ø8 mm) no se malla: perturbación pequeña; su acoplamiento se calcula con el flujo del modo
# propio (fem/rf_acoplos.py). dee_bridge_ entra por el prefijo "dee_"


def dummy_dees(P):
    """Dummy dees (a tierra) en los bordes de las colinas, frente a los bordes del dee: forman los 4 gaps aceleradores
    (diseño IBA de valle profundo, JUAS 2022 p.44–45). Gap constante 'dee.dummy_gap' respecto del borde radial del dee;
    se extienden hasta el borde de la colina (con 2 mm de solape). Solo si dee.dummy_gap > 0."""
    DE = P["dee"]
    g = DE.get("dummy_gap", 0.0)
    if not g:
        return []
    zi = DE["aperture"] / 2
    t = DE.get("dummy_t", DE["plate_t"])
    half_valley = (360.0 / P["pole"]["n_sectors"] - P["pole"]["hill_angle"]) / 2
    r0, r1 = DE["r_min"] + 20, DE["r_max"]
    out = []
    for az in DE["valleys"]:
        for s in (1, -1):
            te = math.radians(az + s * DE["angle"] / 2); th = math.radians(az + s * half_valley)
            ne = (-math.sin(te) * s, math.cos(te) * s); nh = (-math.sin(th) * s, math.cos(th) * s)
            pi = lambda r: (r * math.cos(te) + g * ne[0], r * math.sin(te) + g * ne[1])     # noqa: E731
            po = lambda r: (r * math.cos(th) + 2 * nh[0], r * math.sin(th) + 2 * nh[1])     # noqa: E731
            pts = [pi(r0), pi(r1), po(r1), po(r0)]
            out.append(cq.Workplane("XY").workplane(offset=zi).polyline(pts).close().extrude(t).val())
    return out


def dominio_rf(P, E, log=print):
    parts, info = build(P, E, "render")
    Y, VC, X, DE = P["yoke"], P["vacuum_chamber"], P["exterior"], P["dee"]
    H = Y["height"] / 2
    zA = info["zA"]
    per = 360.0 / P["pole"]["n_sectors"]
    valles = [P["pole"]["hill_center0"] + per / 2 + k * per for k in range(P["pole"]["n_sectors"])]
    dee_az = [round(a) % 360 for a in DE["valleys"]]
    vac = zcyl(VC["inner_radius"], 0, zA)
    for az in valles:
        a = math.radians(az)
        x, y = Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a)
        vac = vac.fuse(zcyl(Y["pump_hole_radius"], 0, H, x, y))
        if round(az) % 360 in dee_az:
            vac = vac.fuse(zcyl(X["rf_housing_d"] / 2 - 20, H, H + X["rf_housing_h"] - 20, x, y))
    vac = vac.clean()
    cond = [p for p in parts if p["name"].startswith(CONDUCTORES) and not p["name"].endswith("_lower")]
    bb = vac.BoundingBox()
    for p in cond:
        b = p["shape"].BoundingBox()
        if b.xmax < bb.xmin or b.xmin > bb.xmax or b.ymax < bb.ymin or b.ymin > bb.ymax or b.zmax < bb.zmin or b.zmin > bb.zmax:
            continue
        vac = vac.cut(p["shape"])
    vac = vac.clean()
    sol = vac.Solids()
    if len(sol) > 1:                                          # quedarse con el volumen conectado principal
        sol = sorted(sol, key=lambda s: -s.Volume())
        log(f"[rf] {len(sol)} sólidos; se toma el mayor ({sol[0].Volume() / 1e9:.4f} m³), "
            f"descartados {[round(s.Volume() / 1e6, 2) for s in sol[1:]]} L")
        vac = sol[0]
    log(f"[rf] dominio de vacío (mitad): {vac.Volume() / 1e9:.4f} m³, válido {vac.isValid()}, "
        f"conductores restados: {len(cond)}")
    return vac, info


def exportar_step(P, E, ruta, log=print, paso_colina=45.0):
    # el polígono fino de la colina con shim (15 mm) multiplica las caras pequeñas y la malla RF no cabe en 7 GB;
    # con 45 mm (más los nodos del shim) el error de forma es < 0.05 mm, irrelevante para el modo de 40 MHz
    import kiube.geometry as G
    G.HILL_STEP_MM = paso_colina
    vac, info = dominio_rf(P, E, log)
    vac_m = vac.scale(0.001)                                  # metros para el solver
    cq.exporters.export(vac_m, str(ruta))
    return ruta, info


if __name__ == "__main__":
    from kiube import load_params
    from kiube.beam import extraction
    P, _ = load_params(); E = extraction(P)
    out = ROOT / "out/fem_rf"; out.mkdir(parents=True, exist_ok=True)
    exportar_step(P, E, out / "dominio_rf_m.step")
