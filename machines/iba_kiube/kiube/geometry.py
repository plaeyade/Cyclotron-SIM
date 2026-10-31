"""Geometría paramétrica CadQuery del Cyclone KIUBE (v0.2).

Arquitectura magnética según TUP04 Fig.1:
  A tapa (return yoke)  ·  B yugo lateral  ·  C sector (disco en 4 piezas de 90°)
  D polo (2 escalones)  ·  E pole insert   ·  F tapón central
Origen: centro de la máquina, plano medio z = 0, mm.
"""
import math
import cadquery as cq
from cadquery import Vector as V

LODS = {"print", "fem_struct", "neutronics", "render"}
NOPRINT = {"neutronics"}                 # láminas µm, fluidos (FEM: solo si se pide)
COSMETIC = {"print", "render"}           # exterior de apariencia (ficha técnica): faldón, tuberías, marcos, soportes
NEUT = {"neutronics"}
STRUCT = {"print", "fem_struct", "neutronics", "render"}


def rz(shape, deg):
    return shape.rotate(V(0, 0, 0), V(0, 0, 1), deg)


def box(cx, cy, cz, lx, ly, lz, rot=0.0):
    b = cq.Solid.makeBox(lx, ly, lz, pnt=V(-lx / 2, -ly / 2, -lz / 2))
    return rz(b, rot).translate(V(cx, cy, cz))


def axis_box(org, ang_deg, s0, s1, w, h):
    """Caja alineada con el eje de un haz en el plano medio: de s0 a s1 a lo largo de la dirección ang_deg desde org,
    ancho w (horizontal, perpendicular al haz) y alto h centrado en z = 0."""
    a = math.radians(ang_deg); L = s1 - s0; sm = (s0 + s1) / 2
    return box(org[0] + sm * math.cos(a), org[1] + sm * math.sin(a), 0, L, w, h, ang_deg)


def zcyl(r, z0, z1, x=0.0, y=0.0):
    return cq.Solid.makeCylinder(r, z1 - z0, pnt=V(x, y, z0), dir=V(0, 0, 1))


def wedge(r_out, a0, a1, z0, z1, r_in=0.0):
    am = math.radians((a0 + a1) / 2)
    A0, A1 = math.radians(a0), math.radians(a1)
    P = lambda r, a: (r * math.cos(a), r * math.sin(a))
    wp = cq.Workplane("XY").workplane(offset=z0)
    if r_in <= 0:
        wp = wp.moveTo(0, 0).lineTo(*P(r_out, A0))
    else:
        wp = wp.moveTo(*P(r_in, A0)).lineTo(*P(r_out, A0))
    wp = wp.threePointArc(P(r_out, am), P(r_out, A1))
    if r_in <= 0:
        wp = wp.close()
    else:
        wp = wp.lineTo(*P(r_in, A1)).threePointArc(P(r_in, am), P(r_in, A0)).close()
    return wp.extrude(z1 - z0).val()


def axis_cyl(r, s0, s1, origin, ang_deg, z=0.0):
    a = math.radians(ang_deg)
    d = V(math.cos(a), math.sin(a), 0)
    return cq.Solid.makeCylinder(r, s1 - s0, pnt=V(origin[0], origin[1], z) + d * s0, dir=d)


def halfplane_cut(phi_deg, outward_sign, offset, z0, z1):
    """Caja que ocupa el semiplano exterior a la recta radial phi, desplazada 'offset' hacia dentro."""
    L = 3000.0
    b = cq.Solid.makeBox(L, L, z1 - z0, pnt=V(-L / 2, -offset if outward_sign > 0 else offset - L, z0))
    return rz(b, phi_deg)


def mirror_z(shape):
    return shape.mirror("XY")


def puller_mordaza(S, en_deg, s0, L, w, h1, z_top, chaflan=0.0, chaflan_ang=45.0):
    """Mordaza superior del puller (v0.32): perfil en el plano (s, z) con s a lo largo de la dirección de extracción
    e_n (en_deg), desde s0 (cara frontal) hasta s0 + 2L + 6; abertura h1; bisel del borde interior frontal."""
    xe = 2 * L + 6.0
    cz = chaflan * math.tan(math.radians(chaflan_ang)) if chaflan > 0 else 0.0
    pts = [(0.0, z_top), (0.0, h1 / 2 + cz), (chaflan, h1 / 2), (xe, h1 / 2), (xe, z_top)] if chaflan > 0 else \
          [(0.0, z_top), (0.0, h1 / 2), (xe, h1 / 2), (xe, z_top)]
    # perfil en XZ (x = s, z = z) extruido simétrico en y (ancho)
    sol = cq.Workplane("XZ").polyline(pts).close().extrude(w / 2, both=True).val()
    a = math.radians(en_deg)
    return rz(sol, en_deg).translate(V(S[0] + s0 * math.cos(a), S[1] + s0 * math.sin(a), 0))


def fuente_marco(P, k):
    """Centro de la chimenea k y dirección de extracción (v0.32): 40.3° antes del borde de entrada de la punta del dee
    que alimenta (fuente k → dee k), rendija tangente (sentido de giro del H⁻)."""
    IS, DE = P["ion_source"], P["dee"]
    edge = DE["valleys"][k] - DE.get("punta_ang", 50.0) / 2
    th = math.radians(edge - IS.get("ang_antes_borde", 40.3))
    r = IS["chimney_offset"]
    S = (r * math.cos(th), r * math.sin(th))
    return S, math.degrees(th) + 90.0


HILL_STEP_MM = 15.0      # paso radial del polígono de la colina con shim (el dominio RF usa uno más grueso)


def hill_footprint(PO, hc, R, z0, z1):
    hw = PO["hill_angle"] / 2
    shim = PO.get("shim_dalpha")
    if not shim or not any(abs(float(x)) > 1e-9 for x in shim):
        w = wedge(R, hc - hw - 4, hc + hw + 4, z0, z1)
        w = w.cut(halfplane_cut(hc + hw, +1, PO["hill_side_offset"], z0 - 1, z1 + 1))
        w = w.cut(halfplane_cut(hc - hw, -1, PO["hill_side_offset"], z0 - 1, z1 + 1))
        return w
    # v0.33: perfil de isocronía (shimming de los bordes de colina). Cada lado sigue la recta desplazada hill_side_offset
    # hacia dentro más un giro Δα(r) por lado [°], interpolado linealmente entre los radios de control shim_r [mm]
    o = PO["hill_side_offset"]; rc_ = [float(x) for x in PO["shim_r"]]; da = [float(x) for x in shim]
    import numpy as _np
    rr = _np.r_[_np.linspace(30.0, R, 90)]
    def ang(r, sgn):
        base = hw - math.degrees(math.asin(min(o / r, 0.99)))
        return hc + sgn * (base + float(_np.interp(r, rc_, da)))
    P = lambda r, a: (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))
    # v0.34: lados poligonales con vértices en los radios de control y cada ~15 mm (el perfil es lineal a tramos; el
    # término asin(o/r) es suave): reproduce el polígono de 90 vértices del ajuste con menos caras; arco exterior exacto
    rr = _np.unique(_np.r_[_np.arange(30.0, R, HILL_STEP_MM), [x for x in rc_ if 30 < x < R], R])
    s1 = [P(r, ang(r, +1)) for r in rr]
    s2 = [P(r, ang(r, -1)) for r in rr[::-1]]
    a1, a0 = ang(R, +1), ang(R, -1)
    wp = cq.Workplane("XY").workplane(offset=z0).polyline(s1).threePointArc(P(R, 0.5 * (a0 + a1)), s2[0]).polyline(s2, includeCurrent=False)
    return wp.close().extrude(z1 - z0).val()


def target_cavity(TG, ttype):
    """(Ø, profundidad) de la cavidad del blanco: Ø fija y profundidad desde el volumen de la ficha
    (target.volumen_ml: F18-700 3.6 ml, 13N 2.9 ml, Ga68 3.7 ml, 11C 55 ml); sin dato, cavity_depth."""
    V = (TG.get("volumen_ml") or {}).get(ttype)
    d = 30.0 if ttype == "C11" else TG["cavity_d"]          # 11C: cámara de gas Ø30 (EST) para 55 ml
    if not V:
        return d, (80.0 if ttype == "C11" else TG["cavity_depth"])
    return d, V * 1000.0 / (math.pi * (d / 2) ** 2)


def build(P, E, lod="neutronics"):
    Y, PO, CO, VC, DE = P["yoke"], P["pole"], P["coil"], P["vacuum_chamber"], P["dee"]
    IS, EX, TG = P["ion_source"], P["extraction"], P["target"]
    parts = []

    def add(name, shape, mat, lods=STRUCT, role="hardware"):
        parts.append(dict(name=name, shape=shape, mat=mat, lods=set(lods), role=role))

    H = Y["height"] / 2
    zA = H - Y["plate_thickness"]            # cara interior de la tapa A (= 264)
    g = PO["gap"] / 2
    zC = PO["valley_floor_z"]                # cara de C / base de D
    n = PO["n_sectors"]
    per = 360.0 / n
    valleys = [PO["hill_center0"] + per / 2 + k * per for k in range(n)]
    pump_xy = [(Y["pump_hole_r_pos"] * math.cos(math.radians(a)),
                Y["pump_hole_r_pos"] * math.sin(math.radians(a))) for a in valleys]

    # ================= yugo A+B (mitad superior), espejo para la inferior
    s, c = Y["side"], Y["corner_chamfer"]
    hs = s / 2
    pts = [(-hs + c, -hs), (hs - c, -hs), (hs, -hs + c), (hs, hs - c),
           (hs - c, hs), (-hs + c, hs), (-hs, hs - c), (-hs, -hs + c)]
    yoke = cq.Workplane("XY").polyline(pts).close().extrude(H).val()
    ch_e = Y.get("edge_chamfer", 0.0)
    if ch_e > 0:                          # chaflán 45° de la arista exterior de la tapa (ficha técnica)
        w0 = cq.Workplane("XY").workplane(offset=H - ch_e).polyline(pts).close().wires().val()
        w1 = cq.Workplane("XY").workplane(offset=H).polyline(pts).close().wires().val().offset2D(-ch_e, "intersection")[0]
        env = cq.Workplane("XY").polyline(pts).close().extrude(H - ch_e).val().fuse(cq.Solid.makeLoft([w0, w1])).clean()
        yoke = yoke.intersect(env)
    yoke = yoke.cut(zcyl(Y["cavity_radius"], -1, zA))
    for pr in E["ports"]:
        w = Y["port_widths"][pr["sys"]]
        yoke = yoke.cut(box(pr["x"], pr["y"], Y["port_height"] / 4 - 0.5, 900, w,
                            Y["port_height"] / 2 + 1, pr["dir"]))
    for (x, y) in pump_xy:
        yoke = yoke.cut(zcyl(Y["pump_hole_radius"], zA - 1, H + 1, x, y))
    # aberturas (2) de T04F2 = inserción radial de las fuentes por el plano medio [VE24 Fig.1]
    for az in IS.get("axis_az", []):
        a = math.radians(az)
        rc = Y["cavity_radius"] + 300
        yoke = yoke.cut(box(rc * math.cos(a), rc * math.sin(a), Y["is_port_h"] / 4 - 0.5, 700,
                            Y["is_port_w"], Y["is_port_h"] / 2 + 1, az))
    yoke_up, yoke_lo = yoke, mirror_z(yoke)
    add("yoke_AB_upper", yoke_up, "steel_1010", role="iron")
    add("yoke_AB_lower", yoke_lo, "steel_1010", role="iron")

    # ================= F: tapón central (cara a ras del polo, T04F1)
    for half, sgn in (("upper", 1), ("lower", -1)):
        plug = zcyl(PO["central_plug_r"], g + PO.get("plug_recess", 0.0), zA)
        add(f"F_plug_{half}", plug if sgn > 0 else mirror_z(plug), "steel_1010", role="iron")

    # ================= C: sector (disco en 4 piezas valle→valle con medios agujeros de bombeo)
    plug_void = zcyl(PO["central_plug_r"], -H, H)
    holes = [zcyl(Y["pump_hole_radius"], -H, H, x, y) for (x, y) in pump_xy]
    d0, Re = E["edge_d0"], E["edge_R"]
    t2, ch = PO["step2_thickness"], PO["step2_chamfer"]
    for k in range(n):
        hc = PO["hill_center0"] + k * per
        a = math.radians(hc)
        C = wedge(PO["sector_C_radius"], hc - per / 2, hc + per / 2, zC, zA).cut(plug_void)
        for h_ in holes:
            C = C.cut(h_)
        # ---------- D: polo de 2 escalones
        step1 = hill_footprint(PO, hc, PO["step1_radius"], g + t2, zC)
        cx, cy = d0 * math.cos(a), d0 * math.sin(a)
        arc = cq.Solid.makeCone(Re - ch, Re, ch, pnt=V(cx, cy, g), dir=V(0, 0, 1)).fuse(
            zcyl(Re, g + ch, g + t2, cx, cy))
        step2 = hill_footprint(PO, hc, PO["step1_radius"], g, g + t2).intersect(arc)
        D = step1.fuse(step2).clean().cut(plug_void)
        envelope = D
        rin, rout = PO["central_plug_r"] - 5, PO["step1_radius"] + 10   # -5: evita arista cuchillo tangente al tapón
        rcm = (rin + rout) / 2
        groove = box(rcm * math.cos(a), rcm * math.sin(a), g + PO["insert_groove_d"] / 2,
                     rout - rin, PO["insert_groove_w"], PO["insert_groove_d"], hc)
        cuts = []
        for e in [e for e in E["edge_crossings"] if e["pole"] == k]:
            nang = math.degrees(math.atan2(e["y"] - cy, e["x"] - cx))
            na = math.radians(nang)
            if PO.get("gc_forma", "ranura") == "cuna":
                # corrector rebajado en cuña (PAT16, US10278277): ángulo de 90°, abertura L10 en el borde, vértice a H10
                # hacia dentro de la colina, en la fracción superior (lado del gap) de la altura de la colina
                n_ = (math.cos(na), math.sin(na)); t_ = (-math.sin(na), math.cos(na))
                H10, L10 = PO["gc_cut_depth"], PO["gc_cut_width"]
                Pt = lambda u, v: (e["x"] + u * n_[0] + v * t_[0], e["y"] + u * n_[1] + v * t_[1])   # noqa: E731
                ext = 40.0                                           # se prolonga hacia fuera de la colina
                pts_c = [Pt(-H10, 0.0), Pt(ext, L10 / 2 + ext + H10 - L10 / 2), Pt(ext, -(L10 / 2 + ext + H10 - L10 / 2))]
                hz = PO.get("gc_alto_frac", 0.5) * (zC - g)
                cuts.append(cq.Workplane("XY").workplane(offset=g - 0.1).polyline(pts_c).close().extrude(hz + 0.1).val())
            else:
                L = PO["gc_cut_depth"] + 150.0
                off = -PO["gc_cut_depth"] + L / 2
                cuts.append(box(e["x"] + off * math.cos(na), e["y"] + off * math.sin(na),
                                (g + zC) / 2, L, PO["gc_cut_width"], zC - g + 0.2, nang))
        D = D.cut(groove)
        for cu in cuts:
            D = D.cut(cu)
        ins = box(rcm * math.cos(a), rcm * math.sin(a),
                  g + PO["insert_recess"] + (PO["insert_groove_d"] - PO["insert_recess"]) / 2,
                  rout - rin, PO["insert_groove_w"], PO["insert_groove_d"] - PO["insert_recess"], hc)
        ins = ins.intersect(envelope)
        for cu in cuts:
            ins = ins.cut(cu)
        for half, f in (("upper", lambda x: x), ("lower", mirror_z)):
            add(f"C_sector_{k}_{half}", f(C), "steel_1010", role="iron")
            add(f"D_pole_{k}_{half}", f(D), "steel_1010", role="iron")
            add(f"E_insert_{k}_{half}", f(ins), "steel_1010", role="iron")

    # ================= bobinas (apoyadas sobre A)
    coil = zcyl(CO["inner_radius"] + CO["radial_width"], CO["z_bottom"], CO["z_bottom"] + CO["height"]).cut(
        zcyl(CO["inner_radius"], CO["z_bottom"] - 1, CO["z_bottom"] + CO["height"] + 1))
    add("coil_upper", coil, "copper_coil_homog", role="coil")
    add("coil_lower", mirror_z(coil), "copper_coil_homog", role="coil")

    # ================= cámara de vacío (anillo polo–bobina, entre tapas A)
    ri, ro = VC["inner_radius"], VC["inner_radius"] + VC["wall"]
    chamber = zcyl(ro, -zA, zA).cut(zcyl(ri, -zA - 1, zA + 1))

    # ================= strippers
    inst = EX.get("strippers_instalados", list(range(len(E["strippers"]))))   # PROP: 4 de serie
    for i, st in enumerate(E["strippers"]):
        if i not in inst:
            continue
        phi = st["phi"]
        fs = EX["foil_size"]
        # brazo inclinado aguas arriba (arm_tilt) para no cruzar el haz del blanco vecino
        da = math.radians(phi + EX["arm_tilt"])
        dx, dy = math.cos(da), math.sin(da)
        px, py = st["x"] + (fs + 2) * dx, st["y"] + (fs + 2) * dy   # brazo empieza fuera de la lámina
        bq = px * dx + py * dy
        t_out = -bq + math.sqrt(bq * bq - (px * px + py * py - ro * ro))
        armang = math.degrees(da)
        arm = box(px + t_out / 2 * dx, py + t_out / 2 * dy, 0, t_out, EX["arm_w"], EX["arm_h"], armang)
        chamber = chamber.cut(arm)
        # borde interior de la lámina en R_ex − ΔR/2 (TR24); la lámina se extiende hacia afuera
        rs_ = math.hypot(st["x"], st["y"]); ap = math.radians(phi)
        rf_c = rs_ - E["dR_turn_mm"] / 2 + fs / 2
        foil = box(rf_c * math.cos(ap), rf_c * math.sin(ap), 0, fs, EX["stripper_foil_um"] / 1000, fs, phi)
        ex_, ey_ = px + t_out * dx, py + t_out * dy
        eph = math.atan2(ey_, ex_)
        rcar = ro + 1 + 20
        zcar = min(53.0, CO["z_bottom"] - 20.0)                     # bajo la bobina (v0.18: bobina en z 59)
        car = box(rcar * math.cos(eph), rcar * math.sin(eph), zcar, 40, 50, 38, math.degrees(eph))  # alto 38 (EST)
        add(f"stripper_arm_{i}", arm, "stainless_304L")
        add(f"stripper_foil_{i}", foil, "graphite_foil", NOPRINT)
        add(f"stripper_carousel_{i}", car, "stainless_304L")

    # ================= blancos
    types = TG["types"]
    for i, tg in enumerate(E["targets"]):
        ttype = types[i] if i < len(types) else "F18"
        org, ang = (tg["x"], tg["y"]), tg["dir"]
        cr = TG["collimator_od"] / 2
        chamber = chamber.cut(axis_cyl(cr, -200, 0, org, ang))
        # v0.27: blanco en la cara del yugo (VE24 Fig. 3) → tubo de haz (vacío) desde la pared exterior de la cámara
        # hasta el colimador, por la escotadura del puerto, entre las bobinas
        if EX.get("blanco_en_cara", True):
            a_ = math.radians(ang); ox_, oy_ = org
            bq = ox_ * math.cos(a_) + oy_ * math.sin(a_)
            d_cam = bq - math.sqrt(max(bq * bq - (ox_ * ox_ + oy_ * oy_ - (ro - 1) ** 2), 0.0))   # distancia hacia atrás a r = ro
            ti, te = EX.get("tubo_haz_d", [40.0, 50.0])
            s0 = -d_cam - 2.0
            tube = axis_cyl(te / 2, s0, -TG["collimator_len"], org, ang).cut(axis_cyl(ti / 2, s0 - 1, -TG["collimator_len"] + 1, org, ang))
            tube = tube.cut(zcyl(ro, -zA, zA))                        # empieza en la cara exterior de la cámara
            chamber = chamber.cut(axis_cyl(ti / 2, s0 - 30, s0 + 30, org, ang))
            add(f"T{i+1}_tubo_haz", tube, "stainless_304L")
        coll = axis_cyl(cr, -TG["collimator_len"], 0, org, ang).cut(
            axis_cyl(EX["collimator_d"] / 2, -TG["collimator_len"] - 1, 1, org, ang))
        add(f"T{i+1}_collimator", coll, "aluminum_6061")
        if ttype == "none":
            continue
        t1 = TG["vac_window_um"] / 1000
        cav_d, cav_L = target_cavity(TG, ttype)
        if ttype == "C11":
            t2w, wmat, fill = 0.5, "aluminum_6061", "n2_20bar"
        elif ttype == "Ga68":                         # VE24: Nb 75 µm a 13 MeV (energía variable) / 250 µm a 18 MeV
            t2w = 0.075 if tg.get("wedge_deg") else 0.25
            wmat, fill = "niobium", "zn68_nitrato_33"
        elif ttype == "N13":
            t2w, wmat, fill = 0.025, "havar", "water_nat"
        else:
            t2w, wmat, fill = TG["tgt_window_um"] / 1000, "havar", "water_o18"
        he, wr = TG["he_gap"], TG["window_d"] / 2
        s1, s2, s3 = t1, t1 + he, t1 + he + t2w
        blen = max(TG["body_len"], s3 + cav_L + 15)
        cav = axis_cyl(cav_d / 2, s3, s3 + cav_L, org, ang)
        if abs(wr - cav_d / 2) < 1e-6:                # misma Ø ventana/cavidad: un solo taladro (evita una cara
            body = axis_cyl(TG["body_d"] / 2, 0, blen, org, ang).cut(axis_cyl(wr, -1, s3 + cav_L, org, ang))   # plana de ancho 0)
        else:
            body = axis_cyl(TG["body_d"] / 2, 0, blen, org, ang).cut(axis_cyl(wr, -1, s3, org, ang)).cut(cav)
        add(f"T{i+1}_vac_window", axis_cyl(wr, 0, s1, org, ang), "titanium", NOPRINT)
        add(f"T{i+1}_he_gap", axis_cyl(wr, s1, s2, org, ang), "helium", NEUT)
        add(f"T{i+1}_tgt_window", axis_cyl(wr, s2, s3, org, ang), wmat, NOPRINT)
        add(f"T{i+1}_{ttype}_fill", cav, fill, NOPRINT)
        add(f"T{i+1}_body", body, "niobium")
        # v0.27: puerta blindada Nirta (HDPE, 18 cm aguas abajo del blanco; ficha IBA en VUB 2024) en los puertos del proyecto
        if i in TG.get("puertas_blindadas", []) and TG.get("puerta_hdpe_mm"):
            pw, ph = TG.get("puerta_wh", [260.0, 260.0])
            Lp = blen + TG["puerta_hdpe_mm"]
            door = axis_box(org, ang, 2.0, Lp, pw, ph).cut(axis_cyl(TG["body_d"] / 2 + 5, 0, blen + 5, org, ang))
            add(f"T{i+1}_puerta_blindada", door, "hdpe")
        if tg.get("wedge_deg"):                      # brida en cuña de energía variable [VE24 Fig.6]
            w18 = ang - tg["wedge_deg"]
            Rw_ = TG["body_d"] / 2 + 12
            wedge_fl = axis_cyl(Rw_, -24, 0, org, ang).intersect(axis_cyl(Rw_ * 3, -30, 40, org, w18).cut(
                axis_cyl(Rw_ * 3, -100, -14, org, w18)))
            wedge_fl = wedge_fl.cut(axis_cyl(cr + 0.2, -60, 5, org, ang)).cut(zcyl(ro + 0.5, -zA, zA))
            add(f"T{i+1}_ve_wedge_flange", wedge_fl, "stainless_304L")

    # pasamuros de los vástagos de las fuentes (el vástago llena el agujero: sin fuga)
    for az in IS.get("axis_az", []):
        chamber = chamber.cut(axis_cyl(IS["stem_radius"], ri - 5, ro + 5, (0, 0), az))
    add("vacuum_chamber", chamber, "stainless_304L")

    # ================= dees + stems coaxiales en los agujeros de bombeo
    for j, az in enumerate(DE["valleys"]):
        zi = DE["aperture"] / 2
        ztop = max(zi + DE["plate_t"], DE.get("top_z", 0.0))     # caja del dee (fem_rf calibra top_z con la RF)
        plate = wedge(DE["r_max"], az - DE["angle"] / 2, az + DE["angle"] / 2, zi, ztop, r_in=DE["r_min"])
        tw = DE["plate_t"]
        # caja del dee: chapa de cobre de espesor plate_t. Se modela MACIZA (sin cavidad interior cerrada, que rompería
        # las mallas de vacío y Monte Carlo) con cobre homogeneizado: ρ = ρ_Cu · V_chapa / V_caja ('copper_dee_caja')
        mat_dee = "copper_dee_caja" if ztop - zi > 3 * tw else "copper_ofe"
        a = math.radians(az)
        stem = zcyl(DE["stem_radius"], ztop, H + 200,
                    Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a))
        for half, f in (("upper", lambda x: x), ("lower", mirror_z)):
            add(f"dee_{j}_{half}", f(plate), mat_dee)
            add(f"dee_stem_{j}_{half}", f(stem), "copper_ofe")

    # ================= región central: puntas de los dees unidas por un cuello (v0.22, TUD03T diap. 5 = mapa Opera de IBA).
    #   Un solo amplificador y un solo cable alimentan los 2 dees [SRG §3.6] → unidos galvánicamente en el centro (modo en
    #   fase, h = 2). Las puntas son placas finas (punta_t) de bordes radiales a ±punta_ang/2 desde r = punta_r hasta r_min,
    #   donde empieza la caja del dee; el cuello (puente_w) une las dos puntas, desplazado puente_off del eje.
    #   z de zi a zi + punta_t: por encima de las chimeneas (|z| ≤ 11) y por debajo del tapón (z ≥ g + plug_recess)
    wb = DE.get("puente_w", 0.0)
    if wb and len(DE["valleys"]) == 2:
        zi = DE.get("punta_z0", DE["aperture"] / 2)                  # semiabertura en la región central (v0.24)
        tt = DE.get("punta_t", DE["plate_t"])
        rt, ha = DE.get("punta_r", 8.0), DE.get("punta_ang", 50.0) / 2
        a0 = DE["valleys"][0]
        br = None
        z_caja = DE["aperture"] / 2
        for az_ in DE["valleys"]:
            tip = wedge(DE["r_min"], az_ - ha, az_ + ha, zi, zi + tt, r_in=rt)
            if zi + tt < z_caja:                                      # escalón de unión punta → caja del dee (r_min − 5 … r_min)
                # llega 3 mm bajo la caja (contacto de superficie, no de arista: una arista común rompe la malla de RF)
                tip = tip.fuse(wedge(DE["r_min"] + 3, az_ - DE["angle"] / 2, az_ + DE["angle"] / 2, zi + tt - 0.5, z_caja,
                                     r_in=DE["r_min"] - 5))
            br = tip if br is None else br.fuse(tip)
        off = DE.get("puente_off", 0.0)
        ao = math.radians(a0 - 90)                                    # dirección perpendicular al eje de los dees
        neck = box(off * math.cos(ao), off * math.sin(ao), zi + tt / 2, wb, 2 * rt + 4, tt, a0 - 90)
        br = br.fuse(neck).clean()
        # v0.32: mordazas del puller unidas a la punta de cada dee (igual que fem/central_region.py, preset kiube)
        ISp = P["ion_source"]
        if ISp.get("puller_lwh"):
            Lp, wp_, hp = ISp["puller_lwh"]; ch_, cha_ = ISp.get("puller_chaflan", [0.0, 45.0])
            for k in range(len(DE["valleys"])):
                S_, en_ = fuente_marco(P, k)
                jaw = puller_mordaza(S_, en_, ISp["chimenea_r"] + ISp["gap_puller"], Lp, wp_, hp, zi + tt, ch_, cha_)
                br = br.fuse(jaw).clean()
        for half, f in (("upper", lambda x: x), ("lower", mirror_z)):
            add(f"dee_bridge_{half}", f(br), "copper_ofe")

    # ================= acoplador de potencia y sintonizador en las carcasas RF superiores (v0.19, EST; KEY22: coupler y
    #   tuner en penetraciones opuestas). Acoplador = lazo inductivo junto a la tapa (máximo de corriente del stem λ/4);
    #   sintonizador = pistón capacitivo de cobre que baja hacia la caja del dee
    RFp, XX = P.get("rf"), P["exterior"]
    if RFp and len(DE["valleys"]) == 2:
        ztop_d = max(DE["aperture"] / 2 + DE["plate_t"], DE.get("top_z", 0.0))
        zcap = H + XX["rf_housing_h"] - 20                            # cara interior de la tapa de la carcasa
        # acoplador (dee 0)
        a = math.radians(DE["valleys"][0]); hx, hy = Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a)
        al = a + math.radians(RFp["acoplador_az"]); ux, uy = math.cos(al), math.sin(al)
        r1, r2 = RFp["acoplador_rho"]; hl = RFp["acoplador_h"]; rw = 4.0
        legs = [zcyl(rw, zcap - hl, zcap, hx + r * ux, hy + r * uy) for r in (r1, r2)]
        foot = box(hx + (r1 + r2) / 2 * ux, hy + (r1 + r2) / 2 * uy, zcap - hl, r2 - r1 + 2 * rw, 2 * rw, 2 * rw, math.degrees(al))
        loop = legs[0].fuse(legs[1]).fuse(foot).clean()
        di, do = RFp["coax_feed_d"]
        fx, fy = hx + r2 * ux, hy + r2 * uy
        # pasamuros + conector dentro de la envolvente del SRG (≤ cima del cáncamo, 730 mm sobre el plano medio)
        feed = zcyl(do / 2 + 3, zcap + 20, zcap + 80, fx, fy).cut(zcyl(do / 2, zcap + 19, zcap + 81, fx, fy))
        feed = feed.fuse(zcyl(di / 2, zcap + 20, zcap + 80, fx, fy)).fuse(
            zcyl(do / 2 + 3, zcap + 80, zcap + 95, fx, fy)).clean()             # conector cerrado arriba (1 sólido)
        add("rf_coupler_loop", loop, "copper_ofe")
        add("rf_coupler_feed", feed, "stainless_304L")
        # sintonizador (dee 1)
        a = math.radians(DE["valleys"][1]); hx, hy = Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a)
        tx, ty = hx + RFp["sintonizador_off"] * math.cos(a), hy + RFp["sintonizador_off"] * math.sin(a)
        zt0 = ztop_d + RFp["sintonizador_gap"]
        if zt0 < zcap - 10:                                           # sintonizador_gap enorme = sin pistón
            add("rf_tuner_plunger", zcyl(RFp["sintonizador_d"] / 2, zt0, zcap, tx, ty), "copper_ofe")
        add("rf_tuner_actuator", zcyl(RFp["sintonizador_d"] / 2 + 15, zcap + 20, zcap + 95, tx, ty), "stainless_304L")

    # ================= dummy dees a tierra en los bordes de las colinas (diseño IBA de valle profundo, JUAS22 p.44–45):
    #   placas de cobre frente a cada borde del dee, a 'dummy_gap' de él (los 4 gaps aceleradores), recortadas contra las
    #   colinas para quedar apoyadas (y a tierra) en ellas
    gd = DE.get("dummy_gap", 0.0)
    if gd:
        zi = DE["aperture"] / 2; td = DE.get("dummy_t", DE["plate_t"])
        hv = (360.0 / PO["n_sectors"] - PO["hill_angle"]) / 2
        hills = [p["shape"] for p in parts if p["name"].startswith(("D_pole_", "C_sector_", "E_insert_"))
                 and p["name"].endswith("_upper")]
        r0, r1 = DE["r_min"] + 20, DE["r_max"]
        for j, az in enumerate(DE["valleys"]):
            for k, s in enumerate((1, -1)):
                te = math.radians(az + s * DE["angle"] / 2); th = math.radians(az + s * (hv + 8.0))
                ne = (-math.sin(te) * s, math.cos(te) * s)
                pdd = [(r * math.cos(te) + gd * ne[0], r * math.sin(te) + gd * ne[1]) for r in (r0, r1)] + \
                      [(r * math.cos(th), r * math.sin(th)) for r in (r1, r0)]
                dd = cq.Workplane("XY").workplane(offset=zi).polyline(pdd).close().extrude(td).val()
                for hsh in hills:
                    if dd.BoundingBox().zmax > hsh.BoundingBox().zmin:
                        dd = dd.cut(hsh)
                dd = dd.clean()
                if len(dd.Solids()) > 1:     # el corte con la colina ensanchada (shim) deja astillas sueltas
                    dd = max(dd.Solids(), key=lambda x: x.Volume())
                for half, f in (("upper", lambda x: x), ("lower", mirror_z)):
                    add(f"dummy_dee_{j}_{k}_{half}", f(dd), "copper_ofe")

    # ================= fuentes PIG (2), inserción RADIAL por el plano medio [VE24 Fig.1, KEY22]
    #   chimenea vertical en el gap del tapón + vástago delgado (gap de colina ±12 mm) + tubo portador
    #   por el valle, la pared de la cámara, entre bobinas y por la abertura (2) del yugo
    for k, az in enumerate(IS["axis_az"]):
        a = math.radians(az)
        if IS.get("chimenea_r"):
            # v0.32: punta de la chimenea (Ø 2·chimenea_r, |z| ≤ 11) en la posición de la región central y vástago recto
            # desde ella hasta el puerto; la rendija (0.6 × 6 mm) mira en la dirección de extracción e_n
            (cx, cy), en_ = fuente_marco(P, k)
            px, py = IS["inner_len_r"] * math.cos(a), IS["inner_len_r"] * math.sin(a)
            beta = math.degrees(math.atan2(py - cy, px - cx)); Ls = math.hypot(px - cx, py - cy)
            rch = IS["chimenea_r"]; sw, sh = IS.get("rendija_wh", [0.6, 6.0])
            tip = zcyl(rch, -g + 1, g - 1, cx, cy)
            e = math.radians(en_)
            tip = tip.cut(box(cx + rch * math.cos(e), cy + rch * math.sin(e), 0, 1.6, sw, sh, en_))
            body = tip.fuse(axis_cyl(IS["body_radius"] * 0.8, 0.5 * rch, Ls, (cx, cy), beta)).clean()   # arranca dentro de la punta (sin contacto tangente)
        else:
            cx, cy = IS["chimney_offset"] * math.cos(a), IS["chimney_offset"] * math.sin(a)
            body = zcyl(IS["body_radius"], -g + 1, g - 1, cx, cy).fuse(
                axis_cyl(IS["body_radius"] * 0.8, IS["chimney_offset"], IS["inner_len_r"], (0, 0), az)).clean()
        r_face = (Y["side"] / 2) / max(abs(math.cos(a)), abs(math.sin(a)))
        stem = axis_cyl(IS["stem_radius"], IS["inner_len_r"], r_face + 60, (0, 0), az).intersect(
            box(0, 0, 0, Y["side"], Y["side"], 4 * IS["stem_radius"]))          # termina a ras de la cara
        if IS.get("chimenea_r"):
            body = body.cut(stem)                                   # el vástago fino ya no es radial: sin solape con el tubo
        add(f"ion_source_{k}", body, "copper_ofe")
        add(f"ion_source_stem_{k}", stem, "stainless_304L")

    # ================= exterior según ficha técnica (DS) y renders oficiales (DSR)
    floor = H - Y["floor_to_top"]
    X = P.get("exterior", {})
    build_exterior(P, E, add, H=H, floor=floor, pump_xy=pump_xy, valleys=valleys, pts=pts)
    from .auxiliares import build_auxiliares          # servicios de la foto real (VE24 p.7): cosmético
    build_auxiliares(P, add, H=H, floor=floor, pts=pts)

    # v0.27: los blancos en la cara del yugo (puertas, colimadores, tubos) mandan sobre las piezas cosméticas de la cara
    funcional = [p_["shape"] for p_ in parts if p_["name"].startswith("T") and
                 any(k in p_["name"] for k in ("_puerta_blindada", "_collimator", "_tubo_haz", "_body"))]
    for p_ in parts:
        if p_["name"].startswith(("port_frame_", "ion_source_cover_", "ion_source_actuator_")):
            sh = p_["shape"]
            bb = sh.BoundingBox()
            for f_ in funcional:
                b2 = f_.BoundingBox()
                if b2.xmax < bb.xmin or b2.xmin > bb.xmax or b2.ymax < bb.ymin or b2.ymin > bb.ymax:
                    continue
                sh = sh.cut(f_)
            p_["shape"] = sh.clean()
    return [p for p in parts if lod in p["lods"]], dict(floor_z=floor, zA=zA, zC=zC)


def build_fem_mag(P, E, sector=True):
    """Modelo magnetostático: hierro + bobinas + dominio de aire conforme.
    Mitad superior (z ≥ 0; simetría de plano medio: B normal en z=0).
    sector=True → 90° de colina a colina (45°–135°) con condición periódica en los cortes.
    Se omiten piezas no magnéticas (Cu, inox, Al, Nb ≡ aire) y las aberturas de bobina."""
    P2 = {k: dict(v) for k, v in P.items()}
    P2["yoke"]["coil_port_az"] = []
    parts, _ = build(P2, E, "fem_struct")
    keep = [p for p in parts if p["role"] in ("iron", "coil") and p["name"].endswith("upper")]
    H = P["yoke"]["height"] / 2
    R_air = 1.6 * P["yoke"]["side"]
    air = zcyl(R_air, 0, H + 0.8 * P["yoke"]["side"])
    if sector:
        S = wedge(R_air + 10, 45.0, 135.0, -1, H + P["yoke"]["side"])
        keep = [dict(p, shape=p["shape"].intersect(S)) for p in keep]
        keep = [p for p in keep if p["shape"].Volume() > 1e-3]
        air = air.intersect(S)
    for p in keep:                       # resta secuencial (el corte contra un compuesto falla con caras compartidas)
        air = air.cut(p["shape"])
    doms = keep + [dict(name="air", shape=air, mat="air", role="air")]
    return doms


def build_exterior(P, E, add, H, floor, pump_xy, valleys, pts):
    """Exterior del Cyclone KIUBE completado con la ficha técnica V6 y sus renders.
    Cotas DSR = proporciones leídas de las renders (±10 %); se declaran en params.yaml."""
    Y, X, DE, TG = P["yoke"], P["exterior"], P["dee"], P["target"]
    s, c = Y["side"], Y["corner_chamfer"]
    hs = s / 2
    rpos = Y["pump_hole_r_pos"]
    # ---- carcasas RF (huecas, rodean el stem) sobre/bajo los valles de los dees; bridas ciegas en los otros
    rh, hh = X["rf_housing_d"] / 2, X["rf_housing_h"]
    dee_az = [round(a) % 360 for a in DE["valleys"]]
    for j, az in enumerate(valleys):
        a = math.radians(az)
        x, y = rpos * math.cos(a), rpos * math.sin(a)
        if round(az) % 360 in dee_az:
            tube = zcyl(rh, H, H + hh - 20, x, y).cut(zcyl(rh - 20, H - 1, H + hh, x, y))
            fl = zcyl(rh + 30, H, H + 20, x, y).cut(zcyl(rh - 20, H - 1, H + 21, x, y))
            cap = zcyl(rh + 10, H + hh - 20, H + hh, x, y)
            body = tube.fuse(fl).fuse(cap).clean()
            add(f"rf_housing_{j}_upper", body, "stainless_304L")
            add(f"rf_housing_{j}_lower", mirror_z(body), "stainless_304L")
        else:
            add(f"blank_flange_{j}_upper", zcyl(rh + 30, H, H + 20, x, y), "stainless_304L")
            pr = X["pump_d"] / 2
            zb = -H - X["pump_h"]
            pump = zcyl(pr + 30, -H - 25, -H, x, y).fuse(zcyl(pr, zb, -H - 25, x, y)).fuse(
                zcyl(pr * 0.55, zb - 60, zb, x, y)).clean()
            pump = pump.cut(zcyl(pr - 8, zb + 15, -H + 1, x, y))          # cuerpo hueco: garganta abierta al vacío
            add(f"diffusion_pump_{j}", pump, "stainless_304L")
    # ---- fuentes de iones por fuera: pieza blanca vertical + actuador horizontal (render oficial, DSR)
    IS = P["ion_source"]
    cw, chh = IS["cover_wh"]
    for k, az in enumerate(IS["axis_az"]):
        a = math.radians(az)
        r_face = hs / max(abs(math.cos(a)), abs(math.sin(a)))
        fx, fy = r_face * math.cos(a), r_face * math.sin(a)
        nrm = (1, 0) if abs(fx) >= abs(fy) else (0, 1)
        sg = 1 if (fx if nrm[0] else fy) > 0 else -1
        nx, ny = sg * nrm[0], sg * nrm[1]
        nang = math.degrees(math.atan2(ny, nx))
        add(f"ion_source_cover_{k}", box(fx + 20 * nx, fy + 20 * ny, 0, 40, cw, chh, nang), "paint_white", COSMETIC)
        act = axis_cyl(IS["actuator_d"] / 2, 40, 40 + IS["actuator_len"], (fx, fy), nang)
        act = act.fuse(axis_cyl(IS["actuator_d"] / 2 + 12, 40, 52, (fx, fy), nang)).clean()
        add(f"ion_source_actuator_{k}", act, "steel_dark", COSMETIC)
    # ---- patas de sección cuadrada con placas
    L, d = X["leg_square"], X["leg_xy"]
    for j, (sx, sy) in enumerate(((1, 1), (-1, 1), (-1, -1), (1, -1))):
        x, y = d * sx, d * sy
        leg = box(x, y, (floor + 25 - H - 20) / 2, L, L, -H - 20 - floor - 25)
        Lp = X.get("leg_plate", L + 120)                     # pie 350 × 350 (SRG-0315)
        leg = leg.fuse(box(x, y, floor + 12.5, Lp, Lp, 25)).fuse(box(x, y, -H - 10, L + 60, L + 60, 20)).clean()
        add(f"leg_{j}", leg, "steel_1010")
    # ---- gatos de elevación: principal + auxiliar + ojal, en las 2 esquinas de elevación
    cmid = (s - c) / math.sqrt(2) + 76
    for j, az in enumerate(Y["lifting_corners"]):
        a = math.radians(az)
        cx, cy = cmid * math.cos(a), cmid * math.sin(a)
        tx, ty = -math.sin(a), math.cos(a)
        rj = X["jack_d"] / 2
        add(f"lifting_jack_{j}", zcyl(rj, -H - 130, H + 110, cx, cy), "steel_1010")
        eye = box(cx, cy, H + 170, 50, 90, 120, az).cut(
            cq.Solid.makeCylinder(22, 60, pnt=cq.Vector(cx, cy, H + 185) - cq.Vector(math.cos(a), math.sin(a), 0) * 30,
                                  dir=cq.Vector(math.cos(a), math.sin(a), 0)))
        add(f"lifting_eye_{j}", eye.translate(cq.Vector(0, 0, 0)), "steel_1010")
        for k, sg in enumerate((1, -1)):
            ax_, ay_ = cx + sg * 140 * tx + 25 * math.cos(a), cy + sg * 140 * ty + 25 * math.sin(a)
            add(f"lifting_aux_{j}_{k}", zcyl(42, -H - 60, H - 120, ax_, ay_), "steel_1010")
    # ---- cáncamo central M48 (SRG): único punto de izado superior; su cima marca la altura total (1850 mm)
    eb = X.get("eyebolt")
    if eb:
        z_top = floor + eb["top_floor"]
        ri_, t_ = eb["d_int"] / 2, eb["t"]
        Rr = ri_ + t_ / 2                                   # radio medio del anillo
        zc_ = z_top - t_ / 2 - Rr
        ring = cq.Solid.makeTorus(Rr, t_ / 2, pnt=cq.Vector(0, 0, zc_), dir=cq.Vector(0, 1, 0))
        z_base = zc_ - Rr - t_ / 2 + 8                       # el anillo apoya en el collar
        collar = zcyl(55, z_base - 40, z_base, 0, 0)
        shank = zcyl(24, H, z_base - 40, 0, 0)               # vástago M48 hasta la tapa
        add("lifting_eyebolt_top", ring.fuse(collar).fuse(shank).clean(), "steel_1010", COSMETIC)
    # ---- faldón de lamas negro bajo el yugo (cosmético)
    inset = 40.0
    poly = [(x * (hs - inset) / hs, y * (hs - inset) / hs) for x, y in pts] + [((pts[0][0]) * (hs - inset) / hs,
                                                                                 pts[0][1] * (hs - inset) / hs)]
    sh = X["skirt_h"]
    for e in range(8):
        (x0, y0), (x1, y1) = poly[e], poly[e + 1]
        Ls = math.hypot(x1 - x0, y1 - y0)
        ang = math.degrees(math.atan2(y1 - y0, x1 - x0))
        n = max(2, int(Ls // 70))
        for q in range(n):
            t = (q + 0.5) / n
            add(f"skirt_{e}_{q}", box(x0 + t * (x1 - x0), y0 + t * (y1 - y0), -H - (sh - 20) / 2, 22, 30, sh - 20, ang),
                "polymer_black", COSMETIC)
        add(f"skirt_rail_{e}", box((x0 + x1) / 2, (y0 + y1) / 2, -H - sh + 10, Ls - 40, 34, 20, ang), "polymer_black", COSMETIC)
    # ---- anillo de tuberías de agua (cosmético): tramos rectos + codos esféricos
    po, rp_ = 70.0, X["pipe_d"] / 2
    zr = -H - sh - 40
    ring = [(x * (hs + po) / hs, y * (hs + po) / hs) for x, y in pts]
    pipe = None
    for e in range(8):
        (x0, y0), (x1, y1) = ring[e], ring[(e + 1) % 8]
        seg = cq.Solid.makeCylinder(rp_, math.hypot(x1 - x0, y1 - y0), pnt=cq.Vector(x0, y0, zr),
                                    dir=cq.Vector(x1 - x0, y1 - y0, 0))
        elb = cq.Solid.makeSphere(rp_, pnt=cq.Vector(x0, y0, zr), angleDegrees1=-90, angleDegrees2=90)
        pipe = seg.fuse(elb) if pipe is None else pipe.fuse(seg).fuse(elb)
    add("water_pipe_ring", pipe.clean(), "copper_ofe", COSMETIC)
    # ---- marcos de las ventanas de extracción y soportes de blancos (cosmético)
    faces = {}
    for i, pr in enumerate(E["ports"]):
        az = math.degrees(math.atan2(pr["y"], pr["x"])) % 360
        k = int(((az + 45) % 360) // 90)               # cara: 0=+x, 1=+y, 2=−x, 3=−y
        faces.setdefault(k, []).append((i, pr))
    for k, prs in faces.items():
        nrm = (math.cos(math.radians(90 * k)), math.sin(math.radians(90 * k)))
        tan = (-nrm[1], nrm[0])
        us = []
        for i, pr in prs:
            # intersección del eje del puerto con el plano de la cara
            dx, dy = math.cos(math.radians(pr["dir"])), math.sin(math.radians(pr["dir"]))
            den = dx * nrm[0] + dy * nrm[1]
            t = (hs - (pr["x"] * nrm[0] + pr["y"] * nrm[1])) / den
            hx, hy = pr["x"] + t * dx, pr["y"] + t * dy
            us.append((hx * tan[0] + hy * tan[1], i, pr, hx, hy))
        u0, u1 = min(u for u, *_ in us) - 120, max(u for u, *_ in us) + 120
        um = (u0 + u1) / 2
        cxf, cyf = hs * nrm[0] + um * tan[0] + 6 * nrm[0], hs * nrm[1] + um * tan[1] + 6 * nrm[1]
        frame = box(cxf, cyf, 0, 12, u1 - u0, 260, math.degrees(math.atan2(nrm[1], nrm[0])))
        for u, i, pr, hx, hy in us:
            w = Y["port_widths"][pr["sys"]]
            frame = frame.cut(box(pr["x"], pr["y"], 0, 900, w, Y["port_height"], pr["dir"]))
        for az_is in P["ion_source"]["axis_az"]:           # la pieza de la fuente tapa ese tramo del marco
            ai = math.radians(az_is)
            if int(((az_is + 45) % 360) // 90) == k:
                rf_ = hs / max(abs(math.cos(ai)), abs(math.sin(ai)))
                cw_, ch_ = P["ion_source"]["cover_wh"]
                frame = frame.cut(box(rf_ * math.cos(ai), rf_ * math.sin(ai), 0, 120, cw_ + 4, ch_ + 4,
                                      math.degrees(math.atan2(nrm[1], nrm[0]))))
        add(f"port_frame_{k}", frame, "steel_dark", COSMETIC)
    # ---- autoblindaje opcional (ficha: piezas 4.645 × 2.770 × 2.575 m) — cavidad interior EST
    if X.get("self_shield"):
        lx, ly, lz = X["self_shield_piece"]
        tw = X.get("self_shield_wall", 410.0)               # CALC: 48.1 t por mitad (SRG) con ρ ≈ 3.2 g/cm³
        cav = box(0, 0, floor + (lz - tw) / 2, lx - 2 * tw, 2 * ly - 2 * tw, lz - tw)
        for k, sg in enumerate((1, -1)):
            half = box(0, sg * ly / 2, floor + lz / 2, lx, ly, lz).cut(cav)
            add(f"self_shield_{k}", half, "selfshield_magnetita", {"render", "neutronics"})
