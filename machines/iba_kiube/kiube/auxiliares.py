"""Equipamiento de servicio visible en la foto de la máquina real (ficha VE24 p.7), ausente en el render oficial:
colector de agua de latón con válvulas de bola, mangueras azules a las bobinas, transmisores/caudalímetros azules,
acometida negra, caja eléctrica, manguera neumática en espiral y cadena portacables sobre la tapa.

Todo es COSMÉTICO (LOD print y render): no entra en neutrónica, FEM ni vacío. Cotas en params.yaml → 'auxiliares',
etiqueta F7 = medidas en la foto con la altura del yugo como escala (±20 %). La foto corresponde a la unidad de
formación de IBA; la disposición real de servicios depende de cada instalación.

Se construye en el marco local de la cara de servicio (n = +x hacia fuera, t = +y a lo largo de la cara) y se gira
a 'cara_servicio_az'.
"""
from __future__ import annotations

import math
import cadquery as cq

V = cq.Vector
COSMETIC = {"print", "render"}


def _cil(r, a, b):
    a, b = V(*a), V(*b)
    return cq.Solid.makeCylinder(r, (b - a).Length, pnt=a, dir=(b - a).normalized())


def _caja(c, l, extra_rot=0.0):
    x, y, z = c
    s = cq.Solid.makeBox(*l, pnt=V(x - l[0] / 2, y - l[1] / 2, z - l[2] / 2))
    return s.rotate(V(x, y, z), V(x, y, z + 1), extra_rot) if extra_rot else s


def _tubo(r, puntos, radio_curva):
    """Tubo de radio r por una poligonal 3D con codos de radio 'radio_curva' (barrido de un círculo)."""
    P = [V(*p) for p in puntos]
    edges, cur = [], P[0]
    for i in range(1, len(P) - 1):
        a, b, c = P[i - 1], P[i], P[i + 1]
        d1, d2 = (b - a).normalized(), (c - b).normalized()
        cosang = max(-1.0, min(1.0, d1.dot(d2)))
        th = math.acos(cosang)
        if th < 1e-6:
            continue
        k = radio_curva * math.tan(th / 2)
        p1, p2 = b - d1 * k, b + d2 * k
        bis = (d2 - d1).normalized()
        mid = b + bis * (radio_curva / math.cos(th / 2) - radio_curva)
        if (p1 - cur).Length > 1e-6:
            edges.append(cq.Edge.makeLine(cur, p1))
        edges.append(cq.Edge.makeThreePointArc(p1, mid, p2))
        cur = p2
    edges.append(cq.Edge.makeLine(cur, P[-1]))
    path = cq.Wire.assembleEdges(edges)
    d0 = (P[1] - P[0]).normalized()
    prof = cq.Wire.makeCircle(r, P[0], d0)
    return cq.Solid.sweep(prof, [], path, True, False, None, "right")


def build_auxiliares(P, add, H, floor, pts):
    A = P.get("auxiliares")
    if not A or not A.get("activo", False):
        return
    Y, X = P["yoke"], P["exterior"]
    hs, c = Y["side"] / 2, Y["corner_chamfer"]
    rot = A["cara_servicio_az"]
    R = lambda s: s.rotate(V(0, 0, 0), V(0, 0, 1), rot) if rot else s          # noqa: E731

    # ---------------------------------------------------------------- colector de agua (latón) + ramales
    rc = A["colector_d"] / 2
    xh = hs + A["colector_sep"]
    zh = A["colector_z"]
    u0, u1 = A["colector_u"]
    assert u1 + rc < hs - c, "el colector invade el chaflán de esquina"
    header = _cil(rc, (xh, u0, zh), (xh, u1, zh))
    header = header.fuse(cq.Solid.makeSphere(rc, pnt=V(xh, u0, zh))).fuse(cq.Solid.makeSphere(rc * 1.15, pnt=V(xh, u1, zh)))
    z_union = zh - 110
    header = header.fuse(_cil(rc, (xh, u1, zh), (xh, u1, z_union)))                     # codo + bajante
    header = header.fuse(_cil(rc * 1.3, (xh, u1, z_union + 30), (xh, u1, z_union)))    # racor
    header = header.clean()
    add("water_header", R(header), "brass_cuzn39", COSMETIC)
    n = int(A["n_ramales"])
    us = [u0 + 40 + i * (u1 - 100 - u0 - 40) / (n - 1) for i in range(n)]
    for k, u in enumerate(((us[0] + us[1]) / 2, (us[-2] + us[-1]) / 2)):          # soportes a la cara, entre ramales
        add(f"water_bracket_{k}", R(_caja(((hs + xh) / 2, u, zh), (xh - hs, 30, 40)).cut(header)), "steel_dark", COSMETIC)
    # acometida negra (EPDM): del racor hasta el suelo y hacia fuera
    rn = A["manguera_negra_d"] / 2
    zf = floor + rn + 5
    add("hose_black_0", R(_tubo(rn, [(xh, u1, z_union), (xh, u1, zf), (xh + 350, u1, zf)], 3 * rn)),
        "hose_epdm", COSMETIC)

    # ramales con válvula de bola, manguera azul y transmisor en los pares
    rr, rm = A["ramal_d"] / 2, A["manguera_d"] / 2
    zl = A["manguera_z_bajo"]
    tw, td, th = A["transmisor_wdh"]
    x_in = P["exterior"]["leg_xy"] + P["exterior"]["leg_square"] / 2 + 40              # libra las patas
    for i, u in enumerate(us):
        zv = zh - rc - 70                                                            # centro de la válvula
        br = _cil(rr, (xh, u, zh), (xh, u, zv - 45))
        br = br.fuse(_caja((xh, u, zv), (34, 34, 44))).clean().cut(header)             # cuerpo de la válvula
        add(f"water_branch_{i}", R(br), "brass_cuzn39", COSMETIC)
        lever = _caja((xh + 17 + 30, u, zv + 14), (60, 12, 8)).fuse(_caja((xh + 17 + 3, u, zv + 8), (8, 12, 14)))
        add(f"water_valve_handle_{i}", R(lever.clean().cut(br)), "polymer_red", COSMETIC)
        hose = _tubo(rm, [(xh, u, zv - 45), (xh, u, zl), (x_in, u, zl)], 3.5 * rm)
        if i % 2 == 0:
            zt = (zv - 45 + zl) / 2 - 40
            tx = _caja((xh, u, zt), (td, tw, th))
            hose = hose.cut(tx)
            add(f"flow_tx_{i}", R(tx), "transmitter_equiv", COSMETIC)
        add(f"hose_blue_{i}", R(hose), "hose_pu", COSMETIC)

    # ---------------------------------------------------------------- caja eléctrica
    cw, ch, cd = A["caja_whd"]
    uc, zt = A["caja_u"], A["caja_z_sup"]
    x0 = hs + A["colector_sep"] - 15                                                 # por fuera del anillo de tuberías
    assert abs(uc) < hs - c, "el soporte de la caja eléctrica debe apoyar en la cara plana"
    caja = _caja((x0 + cd / 2, uc, zt - ch / 2), (cd, cw, ch))
    caja = caja.fuse(_caja((x0 + cd + 4, uc, zt - ch / 2), (8, cw - 40, ch - 40)))      # puerta
    caja = caja.fuse(_caja((x0 + cd + 14, uc + cw / 2 - 45, zt - ch / 2), (12, 18, 70)))  # tirador
    add("elec_box_0", R(caja.clean()), "cabinet_equiv", COSMETIC)
    add("elec_bracket_0", R(_caja(((hs + x0) / 2, uc, zt - 30), (x0 - hs, cw - 60, 24))), "steel_dark", COSMETIC)

    # ---------------------------------------------------------------- manguera neumática en espiral (azul)
    od, dt, Lh = A["espiral_odl"]
    rh = (od - dt) / 2
    ue, xe = A["espiral_u"], hs + A["colector_sep"] + od / 2 + 20
    ze = -H + 60                                                                     # salida en la cara del yugo
    feed = _tubo(dt / 2, [(hs, ue, ze), (xe, ue, ze), (xe, ue, ze - 60), (xe + rh, ue, ze - 100)], 2.5 * dt)
    hel = cq.Wire.makeHelix(A["espiral_paso"], Lh, rh, center=V(xe, ue, ze - 100 - Lh))
    # sección de 12 lados (no círculo): el barrido circular por la hélice da superficies B-spline que se
    # teselan en ~460 k triángulos; el dodecágono da ~70 k con la misma flecha (diferencia invisible a Ø 10 mm)
    c0 = V(xe + rh, ue, ze - 100 - Lh)
    pp = [c0 + V(dt / 2 * math.cos(2 * math.pi * k / 12), 0, dt / 2 * math.sin(2 * math.pi * k / 12)) for k in range(12)]
    prof = cq.Wire.makePolygon(pp + [pp[0]])
    coil = cq.Solid.sweep(prof, [], hel, True, True)
    add("pneumatic_coil", R(feed.fuse(coil).clean()), "hose_pu", COSMETIC)

    # ---------------------------------------------------------------- cadena portacables sobre la tapa
    inset = A["cadena_inset"]
    w, hgt = A["cadena_wh"]
    paso = A["cadena_paso"]
    t_lat = 12.0                                   # ≥ 1.2 mm a 1:10 (FDM: pared mínima 0.8 mm tras la ranura)
    poly = [(x * (hs - inset) / hs, y * (hs - inset) / hs) for x, y in pts]
    segs = []
    for e in range(8):
        (xa, ya), (xb, yb) = poly[e], poly[(e + 1) % 8]
        L = math.hypot(xb - xa, yb - ya)
        ang = math.degrees(math.atan2(yb - ya, xb - xa))
        nx, ny = -(yb - ya) / L, (xb - xa) / L
        mx, my = (xa + xb) / 2, (ya + yb) / 2
        z0 = H
        seg = None
        Le = L + w * math.tan(math.radians(22.5))                                    # inglete en los vértices
        for sg in (1, -1):
            ox, oy = mx + sg * nx * (w / 2 - t_lat / 2), my + sg * ny * (w / 2 - t_lat / 2)
            s = _caja((ox, oy, z0 + hgt / 2), (Le, t_lat, hgt), ang)
            seg = s if seg is None else seg.fuse(s)
        nl = max(2, int(L // paso))
        for q in range(nl):
            t = (q + 0.5) / nl - 0.5
            cx_, cy_ = mx + t * (xb - xa), my + t * (yb - ya)
            seg = seg.fuse(_caja((cx_, cy_, z0 + hgt - 6), (14, w, 12), ang))         # travesaño superior
            seg = seg.fuse(_caja((cx_, cy_, z0 + 5), (14, w, 10), ang))               # travesaño inferior
            for sg in (1, -1):                                                     # junta entre eslabones
                gx = cx_ + 0.5 * (xb - xa) / nl + sg * nx * (w / 2 - 0.5)
                gy = cy_ + 0.5 * (yb - ya) / nl + sg * ny * (w / 2 - 0.5)
                if q < nl - 1:
                    seg = seg.cut(_caja((gx, gy, z0 + hgt / 2 + 12), (3, 3, hgt - 24), ang))
        segs.append(seg.clean())
    for e in range(1, 8):                                  # ingletes: sin solape entre tramos contiguos
        segs[e] = segs[e].cut(segs[e - 1])
    segs[7] = segs[7].cut(segs[0])
    for e, seg in enumerate(segs):
        add(f"cable_chain_{e}", seg.clean(), "polymer_black", COSMETIC)
