#!/usr/bin/env python3
"""Genera el XML de Molflow+ para la cámara de vacío del KIUBE (F67).

Esquema leído del escritor/cargador oficiales (fork iBaer/molflow-phd-fork,
src/IO/WriterXML.cpp y LoaderXML.cpp):
  SimulationEnvironment(type=molflow) > Geometry > Vertices/Facets/Structures
  + MolflowSimuSettings > Gas/TimeSettings/Motion/Parameters.
Unidades internas Molflow: cm; outgassing en Pa·m³/s (= mbar·l/s × 0,1,
confirmado en Interface/FacetAdvParams.cpp: "sumOutgassing * 10.00:
Pa*m3/s -> mbar*l/s").

Geometría (KIUBE-V04, params.yaml):
  cámara r=548 mm, gap 132 mm (efectivo, CALC declarado);
  4 bocas de bombeo r=127 mm @ r=375 mm en la tapa superior (T04F2;
  disposición de bombas en tapa superior: EST declarado);
  fuente H₂ central: disco r=20 mm en z=0, 2 caras (EST declarado).
Desgasificado de pared: 1e-10 mbar·l/s/cm² (EST, FUENTES_VACIO).
Bombeo: sticking s=0,25 en las bocas (EST).
Gas: H₂O 18,015 g/mol (dominante, NAKHAI2022); T=293,15 K.
"""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial import Delaunay

R_CM = 54.8          # radio cámara (cm)
GAP_CM = 13.2        # gap efectivo (cm)
R_PORT_CM = 12.7     # radio boca de bombeo
R_POS_CM = 37.5      # radio posición bocas
N_BOMBAS = 4
S_BOMBA = 0.25
Q_PARED = 1e-10      # mbar·l/s/cm² (EST)
Q_H2_TOTAL = 1e-3    # mbar·l/s (EST)
MASA = 18.015
TEMP = 293.15


def _facet(parent, fid, indices, sticking=0.0, outgassing_pa_m3_s=0.0,
           des_type=1, is2sided=0):
    f = ET.SubElement(parent, "Facet", id=str(fid))
    ET.SubElement(f, "Sticking", constValue=str(sticking), parameterId="-1")
    ET.SubElement(f, "Opacity", constValue="1", parameterId="-1",
                  is2sided=str(is2sided))
    ET.SubElement(f, "Outgassing", constValue=str(outgassing_pa_m3_s),
                  parameterId="-1", desType=str(des_type), desExponent="1",
                  hasOutgassingFile="0", useOutgassingFile="0")
    ET.SubElement(f, "Temperature", value=str(TEMP), accFactor="1")
    ET.SubElement(f, "Reflection", diffusePart="1", specularPart="0",
                  cosineExponent="1", type="0", enableSojournTime="0",
                  sojournFreq="1e13", sojournE="0.99")
    ET.SubElement(f, "Structure", inStructure="0", linksTo="0")
    ET.SubElement(f, "Teleport", target="0")
    ET.SubElement(f, "Motion", isMoving="0")
    rec = ET.SubElement(f, "Recordings")
    ET.SubElement(rec, "Profile", type="0", name="none")
    ET.SubElement(rec, "Texture", hasMesh="0", texDimX="0", texDimY="0",
                  countDes="0", countAbs="0", countRefl="0", countTrans="0",
                  countDir="0", countAC="0")
    ind = ET.SubElement(f, "Indices", nb=str(len(indices)))
    for i, v in enumerate(indices):
        ET.SubElement(ind, "Indice", id=str(i), vertex=str(v))
    return f


def _normal_ok(p0, p1, p2, hacia):
    n = np.cross(p1 - p0, p2 - p0)
    return np.dot(n, hacia) > 0


def generar_xml(path_out: str, cfg: dict | None = None) -> dict:
    """cfg (F84, opcional): {"R_cm", "gap_cm", "bombas": [(x_cm, y_cm, lado ±1, r_cm, s)], "q_pared" (mbar·l/s/cm²),
    "q_fuente" (mbar·l/s), "masa" (g/mol), "temp" (K)}. Sin cfg: cámara KIUBE-V04 (4 bocas en la tapa superior)."""
    global TEMP
    cfg = cfg or {}
    R_CM_, GAP_CM_ = cfg.get("R_cm", R_CM), cfg.get("gap_cm", GAP_CM)
    q_pared, q_fuente, masa = cfg.get("q_pared", Q_PARED), cfg.get("q_fuente", Q_H2_TOTAL), cfg.get("masa", MASA)
    TEMP = cfg.get("temp", TEMP)
    bombas = cfg.get("bombas") or [(R_POS_CM * math.cos(math.pi / 4 + math.pi / 2 * k),
                                    R_POS_CM * math.sin(math.pi / 4 + math.pi / 2 * k), +1, R_PORT_CM, S_BOMBA)
                                   for k in range(N_BOMBAS)]
    return _generar(path_out, R_CM_, GAP_CM_, bombas, q_pared, q_fuente, masa)


def _generar(path_out, R_CM, GAP_CM, bombas, Q_PARED, Q_H2_TOTAL, MASA):
    verts: list[tuple] = []
    facetas: list[tuple] = []  # (indices, kwargs)

    def add_v(p):
        verts.append(tuple(p))
        return len(verts) - 1

    zh = GAP_CM / 2
    # ---------- lateral: 2 anillos × 72 segmentos ----------
    NS = 72
    for zz in (-zh, zh):
        for k in range(NS):
            a = 2 * math.pi * k / NS
            add_v((R_CM * math.cos(a), R_CM * math.sin(a), zz))
    for k in range(NS):
        k2 = (k + 1) % NS
        q = [(k, k2, NS + k2), (k, NS + k2, NS + k)]
        for tri in q:
            p = [np.array(verts[i]) for i in tri]
            centro = (p[0] + p[1] + p[2]) / 3
            hacia = -centro / np.linalg.norm(centro[:2].tolist() + [0]) \
                if np.linalg.norm(centro[:2]) > 0 else np.array([0, 0, 1.0])
            if not _normal_ok(p[0], p[1], p[2], hacia):
                tri = (tri[0], tri[2], tri[1])
            facetas.append((list(tri), {}))

    # ---------- tapas: Delaunay con huecos de bombeo en la tapa de cada bomba ----------
    def malla_tapa(z, huecos):
        # puntos: rejilla polar + bordes de huecos densos
        pts = []
        # F84: borde de la tapa = exactamente los NS vértices del anillo lateral (geometría estanca). Hasta v3.9.0 la
        # rejilla polar terminaba en el último múltiplo de 3 cm (54 cm con R = 54,8): quedaba una corona abierta entre la
        # tapa y la pared lateral por la que el gas escapaba de la cámara (presiones de Molflow bajas en F67).
        for k in range(NS):
            a = 2 * math.pi * k / NS
            pts.append((R_CM * math.cos(a), R_CM * math.sin(a)))
        anillos = np.arange(0, R_CM - 1.5, 3.0)
        for r in anillos:
            n = max(1, int(2 * math.pi * max(r, 1) / 3.0))
            for k in range(n):
                a = 2 * math.pi * k / n + r
                pts.append((r * math.cos(a), r * math.sin(a)))
        for cx, cy, rp in huecos:
            for k in range(48):
                a = 2 * math.pi * k / 48
                pts.append((cx + rp * math.cos(a), cy + rp * math.sin(a)))
        pts = np.array(pts)
        pts = pts[np.hypot(pts[:, 0], pts[:, 1]) <= R_CM + 1e-9]
        tri = Delaunay(pts)
        caras = []
        for simplex in tri.simplices:
            p = pts[simplex]
            c = p.mean(axis=0)
            if np.hypot(c[0], c[1]) > R_CM:
                continue
            if any(np.hypot(c[0] - cx, c[1] - cy) < rp for cx, cy, rp in huecos):
                continue
            # descartar triángulos que crucen un hueco (arista larga)
            lados = [np.linalg.norm(p[i] - p[(i + 1) % 3]) for i in range(3)]
            if max(lados) > 6.0:
                continue
            caras.append(simplex)
        base = len(verts)
        for p in pts:
            add_v((p[0], p[1], z))
        hacia = np.array([0.0, 0.0, -1.0 if z > 0 else 1.0])
        for s_ in caras:
            tri_i = [base + i for i in s_]
            p = [np.array(verts[i]) for i in tri_i]
            if not _normal_ok(p[0], p[1], p[2], hacia):
                tri_i = [tri_i[0], tri_i[2], tri_i[1]]
            facetas.append((tri_i, {}))

    for lado in (-1, +1):
        malla_tapa(lado * zh, [(b[0], b[1], b[3]) for b in bombas if b[2] == lado])

    # ---------- bocas de bombeo: discos con sticking s en su tapa ----------
    for cx, cy, lado, rp, sb in bombas:
        zz = lado * zh
        centro = add_v((cx, cy, zz))
        anillo = []
        for k in range(48):
            a = 2 * math.pi * k / 48
            anillo.append(add_v((cx + rp * math.cos(a), cy + rp * math.sin(a), zz)))
        for k in range(48):
            tri = [centro, anillo[k], anillo[(k + 1) % 48]]
            p = [np.array(verts[i]) for i in tri]
            if not _normal_ok(p[0], p[1], p[2], np.array([0.0, 0.0, -float(lado)])):
                tri = [tri[0], tri[2], tri[1]]
            facetas.append((tri, {"sticking": sb}))

    # ---------- fuente H₂ central: disco r=2 cm en z=0, 2 caras ----------
    centro = add_v((0.0, 0.0, 0.0))
    anillo = []
    for k in range(24):
        a = 2 * math.pi * k / 24
        anillo.append(add_v((2.0 * math.cos(a), 2.0 * math.sin(a), 0.0)))
    a_fuente_cm2 = math.pi * 2.0 ** 2
    # F84: en Molflow la desgasificación de una faceta de 2 caras es el TOTAL de la faceta (emite por ambas caras).
    # Hasta v3.9.0 se dividía entre 2 («mitad por cara») y la fuente inyectaba la mitad del H₂ (absorbido = 50 %).
    q_fuente_pa = Q_H2_TOTAL * 0.1
    for k in range(24):
        tri = [centro, anillo[k], anillo[(k + 1) % 24]]
        p = [np.array(verts[i]) for i in tri]
        if not _normal_ok(p[0], p[1], p[2], np.array([0.0, 0.0, 1.0])):
            tri = [tri[0], tri[2], tri[1]]
        facetas.append((tri, {"outgassing_pa_m3_s": q_fuente_pa / 24,
                              "is2sided": 1}))

    # ---------- desgasificado de pared por área ----------
    final = []
    for indices, kw in facetas:
        p = [np.array(verts[i]) for i in indices]
        area_m2 = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) * 1e-4
        kw2 = dict(kw)
        if "outgassing_pa_m3_s" not in kw2:
            area_cm2 = area_m2 * 1e4
            kw2["outgassing_pa_m3_s"] = Q_PARED * area_cm2 * 0.1
        area_cm2_chk = area_m2 * 1e4
        if area_cm2_chk < 1e-8:
            continue  # descartar facetas degeneradas (segfault en molflow)
        final.append((indices, kw2))

    # ---------- XML ----------
    root = ET.Element("SimulationEnvironment", type="molflow",
                      version="2.11.1")
    geo = ET.SubElement(root, "Geometry")
    vg = ET.SubElement(geo, "Vertices", nb=str(len(verts)))
    for i, (x, y, z) in enumerate(verts):
        ET.SubElement(vg, "Vertex", id=str(i), x=str(x), y=str(y), z=str(z))
    fg = ET.SubElement(geo, "Facets", nb=str(len(final)))
    for i, (indices, kw) in enumerate(final):
        _facet(fg, i, indices, **kw)
    st = ET.SubElement(geo, "Structures", nb="1")
    ET.SubElement(st, "Structure", id="0", name="Camara")
    sim = ET.SubElement(root, "MolflowSimuSettings")
    ET.SubElement(sim, "Gas", mass=str(MASA), enableDecay="0",
                  halfLife="100")
    ts = ET.SubElement(sim, "TimeSettings", timeWindow="0.001",
                       useMaxwellDistr="1", calcConstFlow="1")
    ET.SubElement(ts, "UserMoments", nb="0")
    ET.SubElement(sim, "Motion", type="0")
    ET.SubElement(sim, "Parameters")

    ET.indent(ET.ElementTree(root))
    ET.ElementTree(root).write(path_out, xml_declaration=True,
                               encoding="utf-8")
    return {"vertices": len(verts), "facetas": len(final),
            "bombas": len(bombas), "path": path_out}


if __name__ == "__main__":
    import sys
    info = generar_xml(sys.argv[1] if len(sys.argv) > 1
                       else "/tmp/kiube_camara.xml")
    print(info)
