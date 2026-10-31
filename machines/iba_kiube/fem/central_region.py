"""Región central 3D por elementos finitos + órbitas con física real (parametrizado).

Modelo electrostático 3D de la región central (fuente PIG interna, puller, dees con su abertura,
colinas y tapón central a tierra) resuelto con NGSolve (H1 de orden 2) sobre una malla tetraédrica
GRADUADA: muy fina en el gap fuente–puller y en la rendija (≈0.2 mm), ~1 mm en las narices de los
dees y gruesa (≈6 mm) en el resto, de modo que cabe la estructura completa de los dees.
Se resuelve una sola vez el potencial unitario φ₁ (dees = 1, tierra = 0); el campo de RF es
E(x, t) = −∇φ₁ · V_dee · cos(ω t + φ₀)  (cuasi-estático: λ_RF ≈ 7 m ≫ región central).

Tracking: ecuación de Lorentz relativista (Boris) para H⁻ con E de la malla FE (evaluación directa
en el elemento, sin interpolar en grilla) y B(r, φ, z) del imán:
  - preset 'kiube': mapa Bz(r, φ) del modelo 3D no lineal (fem/magnet_ngsolve.py) o, si no existe,
    ⟨B⟩(r) y flutter digitalizados de KEY22 Fig.3;
  - Br, Bφ por expansión de Maxwell de primer orden en z (∇×B = 0, ∇·B = 0).
Fuente: partículas emitidas desde la superficie de plasma (fondo de la rendija, a potencial de la
fuente; electrodo «imaginario» de Smirnov 2021, §5.5) con energía inicial de ~10 eV.

TODA la geometría está en PRESETS (parametrización completa) → modificar y re-optimizar es cambiar
un número. La geometría real de IBA no es pública: los valores marcados EST son de diseño típico.
Uso:
  python fem/central_region.py --preset kiube [--fases -40 40 5] [--vueltas 8] [--h_gap 0.2]
"""
from __future__ import annotations

import argparse, json, math, pathlib, sys, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

C = 299792458.0
E0 = 1.602176634e-19
M_HMINUS_MEV = 939.294           # m(H⁻) c² ≈ m_p + 2 m_e − afinidad

M_DMINUS_MEV = 1876.124        # m(D⁻) c² = m_d + 2 m_e − afinidad (0.75 eV)
ESPECIES = {"H-": M_HMINUS_MEV, "D-": M_DMINUS_MEV}

PRESETS = {
    # KIUBE: RF y voltaje de documentos (T04F8, TUD03); geometría de la región central EST
    "kiube": dict(
        # v0.22: mapa Opera de la región central de IBA (TUD03T diap. 5): puntas de dee de 50° (±25°) hasta r ≈ 8 mm
        # unidas por un cuello; fuentes gemelas a r ≈ 12.5 mm, 40° antes del borde de entrada del dee (aquí: 39.4°
        # = (chimenea + gap + puller)/r). La 2.ª fuente (gemela, a 180°) es de reserva: se simula una activa
        f_rf=40.646e6, V_dee=40.0e3, dee_az=[90.0, 270.0], dee_angle=50.0,
        hill_angle=45.0, hill_gap_half=12.0, plug_r=52.0,
        ap_half=8.0, plate_t=4.0, nose_r=8.0,
        fuentes=[dict(especie="H-", h=2, dee=1, src_r=12.5, src_dtheta=0.0)],
        # slit_h 6 (v0.25, CALC): máximo de corriente útil ∝ altura × transmisión (barrido 3–8 mm,
        # validation/region_central_rendija.json); con 3 mm la transmisión sube a 55 % pero la fuente emite la mitad
        chimney_r=3.0, slit_w=0.6, slit_h=6.0, plasma_depth=2.6,
        # gap_sp 2.8 (v0.23, CALC): óptimo de transmisión del barrido 2.0–3.2 mm (validation/region_central_gap.json,
        # meseta 2.6–3.0); con él la fuente queda a (3 + 2.8 + 3)/12.5 rad = 40.3° del borde, igual que en el mapa de IBA.
        # Campo en el gap 14.3 MV/m (1.7 Kp). El C18/9 usa ≈1.5 mm a 32 kV (JUAS22); escalado a 40 kV (2.0 mm) pierde 7 %
        # puller (v0.24, CALC; validation/region_central_puller.json): abertura vertical 10 mm (barrido 3–12: con < 7 el
        # haz choca con las mordazas; con 12 crecen las pérdidas verticales), ancho 7 (5–10) y profundidad 3 (2–5).
        # ap_half 8 = punta_z0 del CAD: con ±14 la región central no acelera (aceptancia 20–30°)
        gap_sp=2.8, puller_len=3.0, puller_w=7.0, puller_slit_h=10.0,
        # v0.29: cuello entre las puntas (TUD03T diap. 5, igual que el CAD) y giro de la rendija de la chimenea
        puente_w=3.5, puente_off=1.25, slit_rot=0.0,
        # v0.31: forma de las mordazas del puller. Bisel de 5 mm a 45° en el borde interior frontal (CALC,
        # validation/region_central_v31.json): +38 % de corriente a 6 vueltas (plasma plano) y +37 % con el menisco
        # autoconsistente; divergencia vertical a la salida del puller 219 → 110 mrad; campo máximo 37 → 14 MV/m.
        # El escalón de la abertura (puller_h2/puller_d1) no ayuda
        puller_h2=0.0, puller_d1=0.0, puller_chaflan=5.0, puller_chaflan_ang=45.0,
        # v0.29: como el CAD — caja del dee desde r_min = 60 a ±14 (aperture 28) con dee.angle 39°, tapón rebajado 15
        r_caja=60.0, ap_caja=14.0, dee_angle_caja=39.0, plug_recess=15.0,
        dummy=False, dummy_rmin=40.0, dummy_angle=30.0,
        R_dom=140.0, H_dom=32.0, R_fino=40.0,
        h_gap=0.2, h_nose=1.0, h_mid=6.0, h_far=6.0,
        B_source="key22", note="RF/V: T04F8/TUD03; campo: KEY22 (medido); región central: dee 50°, punta r 8 y fuente r 12.5 del mapa Opera de IBA (TUD03T); chimenea, rendija y puller EST"),
    # Cyclone 18/9 — datos publicados: lista de ciclotrones Cyclotrons 2004 (CM06) y 2001 (CU-08):
    #   f_RF 42 MHz (41.8), armónicos 2 (p) / 4 (d), V_dee 32 kV, dee 30°, ⟨B⟩ 1.35 T (valle 0.35, colina 1.9),
    #   gap de colina 30 mm, R_ext 0.48 m, 2 fuentes PIG internas, radio de inyección 0.03 m.
    #   Región central (fuente, puller, narices, dummy dees): EST, tipo IBA (CAS 2017, Figs. 35–37).
    "c18_9": dict(
        f_rf=42.0e6, V_dee=32.0e3, dee_az=[90.0, 270.0], dee_angle=30.0,
        hill_angle=58.0, hill_gap_half=15.0, plug_r=50.0,
        ap_half=10.0, plate_t=4.0, nose_r=7.0,
        # radios: optimizador de centrado (--optimizar); azimut +4°: barrido de ajuste longitudinal (CAS Fig. 37);
        # corte_ext: recorte de la chimenea D⁻ para el paso del H⁻ (CAS Fig. 37), 0 = sin recorte
        fuentes=[dict(especie="H-", h=2, dee=1, src_r=12.8, src_dtheta=4.0),
                 dict(especie="D-", h=4, dee=0, src_r=24.0, src_dtheta=4.0, corte_ext=0.0)],
        chimney_r=3.0, slit_w=0.6, slit_h=6.0, plasma_depth=2.6,
        gap_sp=1.5, puller_len=3.0, puller_w=7.0, puller_slit_h=5.0,
        dummy=True, dummy_rmin=45.0, dummy_angle=30.0,
        R_dom=480.0, H_dom=30.0, R_fino=40.0, mitad=True,
        h_gap=0.2, h_nose=1.0, h_mid=5.0, h_far=18.0,
        B_source="c18_9", B_hill=1.9, B_valley=0.35, B_avg=1.35, bump_G=250.0, bump_sigma=30.0,
        note="CM06/CU-08 (f, h, V, ángulo de dee, campos, gap); gap chimenea–puller 1.5 mm (JUAS22 p.45); resto de la región central EST"),
}


# ----------------------------------------------------------------------------- geometría
def _sector_prism(R, a0, a1, z0, z1, n=48):
    from netgen.occ import WorkPlane, Axes, Pnt, Z, X
    pts = [(0.0, 0.0)] + [(R * math.cos(math.radians(a0 + (a1 - a0) * k / n)),
                           R * math.sin(math.radians(a0 + (a1 - a0) * k / n))) for k in range(n + 1)]
    wp = WorkPlane(Axes((0, 0, z0), n=Z, h=X)).MoveTo(*pts[0])
    for p in pts[1:]:
        wp = wp.LineTo(*p)
    return wp.Close().Face().Extrude(z1 - z0)


def mordazas_puller(pr, z_top):
    """Mordazas del puller en el marco local de la fuente: x a lo largo de e_n con x = 0 en la cara frontal
    (s = chimenea + gap), y a lo ancho (u), z vertical. v0.31: abertura en dos escalones y bisel del borde interior.
      puller_slit_h    abertura vertical en la entrada [mm]
      puller_h2, puller_d1  abertura más adentro (a partir de x = puller_d1); 0 = sin escalón
      puller_chaflan, puller_chaflan_ang  bisel del borde interior frontal (catetos: chaflán en x, chaflán·tan(ang) en z)"""
    from netgen.occ import Box, Pnt, Vec, Axis, Y
    L, pw, h1 = pr["puller_len"], pr["puller_w"], pr["puller_slit_h"]
    h2, d1 = pr.get("puller_h2", 0.0), pr.get("puller_d1", 0.0)
    c, ang = pr.get("puller_chaflan", 0.0), pr.get("puller_chaflan_ang", 45.0)
    xe = 2 * L + 6.0
    out = []
    for sg in (1, -1):
        def caja(x0, x1, hh):
            z0, z1 = (hh / 2, z_top) if sg > 0 else (-z_top, -hh / 2)
            return Box(Pnt(x0, -pw / 2, z0), Pnt(x1, pw / 2, z1))
        if h2 and d1:
            j = caja(0.0, d1, h1) + caja(d1, xe, h2)
        else:
            j = caja(0.0, xe, h1)
        if c > 0:
            cz = c * math.tan(math.radians(ang))
            a = math.hypot(c, cz) * 2
            # prisma triangular: bloque girado que deja un corte con catetos c (en x) y cz (en z) en la arista interior
            th = math.degrees(math.atan2(cz, c))
            zin = sg * h1 / 2
            if sg > 0:   # recta de (0, zin + cz) a (c, zin); se quita lo que queda hacia la arista (abajo-izquierda)
                blk = Box(Pnt(-a, -pw / 2 - 1, -a), Pnt(a, pw / 2 + 1, 0.0))
                cut = blk.Rotate(Axis((0, 0, 0), Y), th).Move(Vec(0, 0, zin + cz))
            else:        # espejo en z
                blk = Box(Pnt(-a, -pw / 2 - 1, 0.0), Pnt(a, pw / 2 + 1, a))
                cut = blk.Rotate(Axis((0, 0, 0), Y), -th).Move(Vec(0, 0, zin - cz))
            j = j - cut
        out.append(j)
    return out


def build(pr):
    """Devuelve (shape, info) con nombres de caras: 'tierra', 'dee', 'exterior'."""
    from netgen.occ import Cylinder, Box, Pnt, Vec, Z, Glue, Axes, Axis, X, OCCGeometry
    R, H = pr["R_dom"], pr["H_dom"]
    # mitad=True: solo z ≥ 0 con simetría de plano medio (E_z = 0 en z = 0, condición natural)
    dom = Cylinder(Pnt(0, 0, 0 if pr.get("mitad") else -H), Z, r=R, h=H if pr.get("mitad") else 2 * H)
    g = pr["hill_gap_half"]
    ground = []
    per = 90.0
    for k in range(4):                                          # colinas (tierra) |z| ≥ g
        hc = 45.0 + k * per
        for s in (1, -1):
            z0, z1 = (g, H + 1) if s > 0 else (-H - 1, -g)
            ground.append(_sector_prism(R + 5, hc - pr["hill_angle"] / 2, hc + pr["hill_angle"] / 2, z0, z1))
    gp = g + pr.get("plug_recess", 0.0)                         # v0.29: tapón rebajado (CAD: 12 + 15 = 27)
    for s in (1, -1):                                           # tapón central
        z0, z1 = (gp, H + 1) if s > 0 else (-H - 1, -gp)
        ground.append(Cylinder(Pnt(0, 0, z0), Z, r=pr["plug_r"], h=z1 - z0))
    # dees: dos placas (superior/inferior) por dee; abertura |z| < ap_half
    dees = []
    ah, t = pr["ap_half"], pr["plate_t"]
    rc_, ac_ = pr.get("r_caja", 0.0), pr.get("ap_caja", ah)
    for az in pr["dee_az"]:
        for s in (1, -1):
            z0, z1 = (ah, ah + t) if s > 0 else (-ah - t, -ah)
            Rt = rc_ if rc_ else R + 5
            w = _sector_prism(Rt, az - pr["dee_angle"] / 2, az + pr["dee_angle"] / 2, z0, z1)
            w = w - Cylinder(Pnt(0, 0, -H - 2), Z, r=pr["nose_r"], h=2 * H + 4)
            dees.append(w)
            if rc_:   # v0.29: caja del dee (r > r_caja) a ±ap_caja con su ángulo, y escalón punta → caja (como el CAD)
                ac = pr.get("dee_angle_caja", pr["dee_angle"])
                zc0, zc1 = (ac_, ac_ + t) if s > 0 else (-ac_ - t, -ac_)
                caja = _sector_prism(R + 5, az - ac / 2, az + ac / 2, zc0, zc1) - Cylinder(Pnt(0, 0, -H - 2), Z, r=rc_, h=2 * H + 4)
                ze0, ze1 = (ah + t - 0.5, ac_ + 0.5) if s > 0 else (-ac_ - 0.5, -ah - t + 0.5)
                esc = (_sector_prism(rc_ + 3, az - ac / 2, az + ac / 2, ze0, ze1)
                       - Cylinder(Pnt(0, 0, -H - 2), Z, r=rc_ - 5, h=2 * H + 4))
                dees += [caja, esc]
    # cuello que une las puntas de los 2 dees en el centro (v0.29; TUD03T diap. 5: 3.5 mm, desplazado 1.25 mm; = CAD)
    if pr.get("puente_w", 0.0) and len(pr["dee_az"]) == 2:
        a0 = math.radians(pr["dee_az"][0]); ao = a0 - math.pi / 2
        L = 2 * pr["nose_r"] + 4.0; wb = pr["puente_w"]; off = pr.get("puente_off", 0.0)
        for s in (1, -1):
            z0, z1 = (ah, ah + t) if s > 0 else (-ah - t, -ah)
            nb = Box(Pnt(-wb / 2, -L / 2, z0), Pnt(wb / 2, L / 2, z1)).Rotate(Axis((0, 0, 0), Z), math.degrees(ao))
            dees.append(nb.Move(Vec(off * math.cos(ao), off * math.sin(ao), 0)))
    # dummy dees: electrodos a tierra en los valles de las fuentes, misma abertura que los dees (EST)
    if pr.get("dummy"):
        for az in pr["dee_az"]:
            for s_ in (1, -1):
                z0, z1 = (ah, ah + t) if s_ > 0 else (-ah - t, -ah)
                w = _sector_prism(R + 5, az + 90 - pr["dummy_angle"] / 2, az + 90 + pr["dummy_angle"] / 2, z0, z1)
                ground.append(w - Cylinder(Pnt(0, 0, -H - 2), Z, r=pr["dummy_rmin"], h=2 * H + 4))
    # fuentes: chimenea a radio src_r ANTES del borde de entrada de su dee (sentido antihorario del ion
    # negativo con Bz > 0), con la rendija en la dirección TANGENTE; puller unido a ese dee
    chims, pl, fu = [], [], []
    for F in pr["fuentes"]:
        azd = pr["dee_az"][F["dee"]]
        edge = azd - pr["dee_angle"] / 2
        rc = pr["chimney_r"]
        reach = rc + pr["gap_sp"] + pr["puller_len"]
        th_s = math.radians(edge) - reach / F["src_r"] + math.radians(F.get("src_dtheta", 0.0))
        S = F["src_r"] * np.array([math.cos(th_s), math.sin(th_s)])
        e_t = np.array([-math.sin(th_s), math.cos(th_s)]); e_r = np.array([math.cos(th_s), math.sin(th_s)])
        rho = math.radians(F.get("slit_rot", pr.get("slit_rot", 0.0)))   # v0.29: giro de la rendija (+ = hacia fuera)
        e_n = math.cos(rho) * e_t + math.sin(rho) * e_r
        pull_front = S + (rc + pr["gap_sp"]) * e_n
        pul_face = pull_front + pr["puller_len"] * e_n
        chim = Cylinder(Pnt(S[0], S[1], -H - 1), Z, r=rc, h=2 * H + 2)
        ang_n = math.degrees(math.atan2(e_n[1], e_n[0]))
        s_cav = pr.get("menisco_s_cav", pr["plasma_depth"])     # v0.30: con menisco, cavidad honda + plasma penalizado
        depth = rc - s_cav + 0.5
        slit = Box(Pnt(0, -pr["slit_w"] / 2, -pr["slit_h"] / 2), Pnt(depth + 1.0, pr["slit_w"] / 2, pr["slit_h"] / 2))
        slit = slit.Rotate(Axis((0, 0, 0), Z), ang_n).Move(Vec(S[0] + s_cav * e_n[0], S[1] + s_cav * e_n[1], 0))
        chim = chim - slit
        if F.get("corte_ext", 0.0) > 0:                       # CAS 2017 Fig. 37: chimenea D⁻ recortada para dar paso al H⁻
            e_r = np.array([math.cos(th_s), math.sin(th_s)])
            c0 = S + (rc - F["corte_ext"]) * e_r
            cut = Box(Pnt(0, -20, -H - 5), Pnt(20, 20, H + 5)).Rotate(Axis((0, 0, 0), Z), math.degrees(th_s)).Move(Vec(c0[0], c0[1], 0))
            chim = chim - cut
        chims.append(chim)
        for b_ in mordazas_puller(pr, ah + t):
            pl.append(b_.Rotate(Axis((0, 0, 0), Z), ang_n).Move(Vec(pull_front[0], pull_front[1], 0)))
        fu.append(dict(especie=F["especie"], h=F["h"], S=S.tolist(), e_n=e_n.tolist(), edge_deg=edge,
                       slit_center=(S + pr["plasma_depth"] * e_n).tolist()))
    for so in ground:
        for f in so.faces:
            f.name = "tierra"
    for so in chims:
        for f in so.faces:
            f.name = "tierra"
    for so in dees + pl:
        for f in so.faces:
            f.name = "dee"
    for f in dom.faces:
        f.name = "exterior"
    if pr.get("mitad"):
        dom.faces.Min(Z).name = "simetria"
    air = dom
    for so in ground + chims + dees + pl:
        air = air - so
    air.mat("vacio")
    info = dict(fuentes=fu, S=fu[0]["S"], e_n=fu[0]["e_n"], edge_deg=fu[0]["edge_deg"], slit_center=fu[0]["slit_center"])
    build.groups = dict(colinas_tapon=ground[:10], dummy_dees=ground[10:], chimeneas=chims, dees=dees, pullers=pl)
    return air, info


def chi_menisco(pr, info):
    """Indicador del plasma: dentro de la rendija y por detrás de s_m(z) = c0 + c1 z² + c2 z⁴ + c3 z⁶ (fem/menisco.py)."""
    from ngsolve import x, y, z, IfPos
    F = info["fuentes"][0]; S = F["S"]; en = F["e_n"]; et = (-en[1], en[0])
    c = pr["_menisco"]
    s_ = (x - S[0]) * en[0] + (y - S[1]) * en[1]; u_ = (x - S[0]) * et[0] + (y - S[1]) * et[1]
    sm = c[0] + c[1] * z ** 2 + c[2] * z ** 4 + c[3] * z ** 6
    if len(c) > 4:
        sm = sm + c[4] * (2 * u_ / pr["slit_w"]) ** 2
    hw, hh = pr["slit_w"] / 2 + 0.05, pr["slit_h"] / 2 + 0.05
    s_cav = pr.get("menisco_s_cav", pr["plasma_depth"])
    return (IfPos(sm - s_, 1.0, 0.0) * IfPos(hw - u_, 1.0, 0.0) * IfPos(hw + u_, 1.0, 0.0) * IfPos(hh - z, 1.0, 0.0) *
            IfPos(hh + z, 1.0, 0.0) * IfPos(s_ - s_cav + 0.3, 1.0, 0.0))


def mesh_params(pr, info):
    """Malla graduada: h_gap en cada gap fuente–puller y rendija; h_nose en las narices de los dees y la
    zona central (r < R_fino); h_mid hasta r ≈ 150 mm; h_far en el resto de la estructura de los dees."""
    from netgen.meshing import MeshingParameters
    mp = MeshingParameters(maxh=pr["h_far"], grading=0.25)
    for F in info["fuentes"]:
        S, en = np.array(F["S"]), np.array(F["e_n"])
        tv = np.array([-en[1], en[0]])
        span = pr["chimney_r"] + pr["gap_sp"] + pr["puller_len"] + 1.0
        for s_ in np.linspace(-0.5, span, 24):
            for tz in np.linspace(-pr["slit_h"] / 2 - 1, pr["slit_h"] / 2 + 1, 7):
                for tt in np.linspace(-pr["puller_w"] / 2, pr["puller_w"] / 2, 5):
                    P_ = S + s_ * en + tt * tv
                    mp.RestrictH(float(P_[0]), float(P_[1]), float(tz), pr["h_gap"])
        if pr.get("menisco_s_cav") is not None:              # v0.30: canal de la rendija resuelto como en fem/menisco.py
            hc = pr.get("h_canal", 0.08)
            for s_ in np.arange(pr["menisco_s_cav"] - 0.1, pr["chimney_r"] + 0.6, hc * 1.2):
                for tt in np.arange(-pr["slit_w"] / 2 - 0.1, pr["slit_w"] / 2 + 0.11, hc * 1.2):
                    for tz in np.arange(-pr["slit_h"] / 2 - 0.3, pr["slit_h"] / 2 + 0.31, hc * 1.2):
                        P_ = S + s_ * en + tt * tv
                        mp.RestrictH(float(P_[0]), float(P_[1]), float(tz), hc)
    Rf = pr.get("R_fino", 36.0)
    for xx in np.arange(-Rf, Rf + 0.1, 3.0):
        for yy in np.arange(-Rf, Rf + 0.1, 3.0):
            if xx * xx + yy * yy < Rf ** 2:
                for zz in (-pr["ap_half"], 0.0, pr["ap_half"]):
                    mp.RestrictH(float(xx), float(yy), float(zz), pr["h_nose"])
    if pr["R_dom"] > 160:
        for rr in np.arange(Rf, 150.1, 8.0):
            for th in np.arange(0, 360, max(2.0, math.degrees(8.0 / rr))):
                c_, s_ = math.cos(math.radians(th)), math.sin(math.radians(th))
                mp.RestrictH(float(rr * c_), float(rr * s_), 0.0, pr["h_mid"])
    return mp


def _v(p):
    from netgen.occ import Vertex
    return Vertex(p)


# ----------------------------------------------------------------------------- campo magnético
class MagField:
    def __init__(self, pr, P):
        self.n = 4
        Mm = pr.get("_mass", M_HMINUS_MEV)
        if pr["B_source"] == "c18_9":
            # ⟨B⟩(r) isócrono para la especie (B0·γ(r)) + joroba central (CAS 2017 §3.3.3: «unos cientos de
            # gauss»; EST) y modulación de 4 colinas: escalón suavizado colina/valle (CM06: 1.9 / 0.35 T) con
            # ancho de borde ≈ medio gap / r → el flutter crece desde 0 en el centro, como en un imán real
            w_rev = 2 * math.pi * pr["f_rf"] / pr["h"]
            B0 = w_rev * Mm * 1e6 / C ** 2
            self.r = np.arange(0.0, pr["R_dom"] + 0.5, 0.5)
            self.ph = np.arange(45.0, 135.01, 0.5)
            g = 1 / np.sqrt(1 - (w_rev * self.r * 1e-3 / C) ** 2)
            Bavg = B0 * g * (1 + pr["bump_G"] * 1e-4 / B0 * np.exp(-self.r ** 2 / (2 * pr["bump_sigma"] ** 2)))
            th = np.radians(self.ph)[None, :]
            wid = np.clip(pr["hill_gap_half"] / np.maximum(self.r, 1.0), 0.02, 3.0)[:, None]
            Sm = np.zeros((len(self.r), len(self.ph)))
            for hc in (45.0, 135.0, -45.0, 225.0):
                a0, a1 = math.radians(hc - pr["hill_angle"] / 2), math.radians(hc + pr["hill_angle"] / 2)
                Sm += 0.5 * (np.tanh((th - a0) / wid) - np.tanh((th - a1) / wid))
            prof = pr["B_valley"] + (pr["B_hill"] - pr["B_valley"]) * np.clip(Sm, 0, 1)
            self.Bz = Bavg[:, None] * prof / prof.mean(axis=1, keepdims=True)
            Fl = ((self.Bz - self.Bz.mean(1, keepdims=True)) ** 2).mean(1) / self.Bz.mean(1) ** 2
            self.flutter = (self.r, Fl)
            self.kind = "map"
            self.src = (f"modelo C18/9: B0 isócrono {B0:.4f} T para {pr.get('_especie', 'H-')} (h = {pr['h']}) + "
                        f"joroba {pr['bump_G']:.0f} G + colinas {pr['B_hill']}/{pr['B_valley']} T suavizadas; "
                        f"flutter a r = 0.3 m: {np.interp(300, self.r, Fl):.3f}")
            return
        if pr["B_source"] == "fem_map":
            cands = sorted((ROOT / "out/fem_mag").glob("ngsolve*/bz_midplane_*.npz"))
            sel = ROOT / "out/fem_mag/campo_calibrado.npz"
            f = sel if sel.exists() else (cands[-1] if cands else None)
            if pr.get("B_npz"):                                  # v0.35: mapa y nivel explícitos (enlace con fem/aceleracion)
                f = pathlib.Path(pr["B_npz"])
            if f is not None:
                d = np.load(f)
                self.r = d["r_m"] * 1e3; self.ph = d["phi_deg"]; Bz = d["Bz_T"]
                # escala para isocronía en el centro: ⟨B⟩(r→0) = B0 del preset (la NI del modelo es EST)
                Bm = np.nanmean(Bz, axis=1)
                B0 = 2 * math.pi * pr["f_rf"] / pr["h"] * Mm * 1e6 * E0 / C ** 2 / E0
                k = B0 / np.nanmean(Bm[(self.r > 20) & (self.r < 80)])
                if pr.get("B_escala"):
                    k = float(pr["B_escala"])
                self.Bz = np.nan_to_num(Bz * k, nan=float(np.nanmean(Bz) * k))
                self.src = f"mapa FE {f.name} × {k:.4f} (⟨B⟩ = B0 isócrono a r = 20–80 mm)"
                self.kind = "map"
                return
        B0 = 2 * math.pi * pr["f_rf"] / pr["h"] * Mm * 1e6 / C ** 2
        self.B0 = B0
        csv = ROOT / "validation/kiube_Bavg_flutter_KEY22.csv"
        if pr["B_source"] in ("fem", "key22") and csv.exists():
            d = np.genfromtxt(csv, delimiter=",", skip_header=4)
            rB = d[~np.isnan(d[:, 1])]; rF = d[~np.isnan(d[:, 2])]
            self.r = np.arange(0, 140.5, 0.5)
            Bm = np.interp(self.r, rB[:, 0] * 1e3, rB[:, 1] * 1e-4)
            Fr = np.interp(self.r, rF[:, 0] * 1e3, rF[:, 2])
            k = B0 / Bm[(self.r > 20) & (self.r < 80)].mean()
            self.ph = np.arange(45.0, 135.01, 1.0)
            amp = np.sqrt(2 * Fr)                                  # modulación sinusoidal: F = a²/2
            self.Bz = (k * Bm)[:, None] * (1 + amp[:, None] * np.cos(np.radians(4 * (self.ph[None, :] - 45.0))))
            self.kind = "map"
            self.src = f"KEY22 Fig.3: ⟨B⟩(r) con joroba central + flutter F(r) (modulación de 4 sectores), × {k:.4f} → B0 isócrono"
            return
        self.kind = "iso"; self.src = f"isócrono uniforme B0 = {B0:.4f} T"

    def __call__(self, X):
        """X: (N,3) en mm → B (N,3) en T."""
        x, y, z = X[:, 0], X[:, 1], X[:, 2]
        r = np.hypot(x, y) + 1e-9
        ph = np.degrees(np.arctan2(y, x)) % 360
        if self.kind == "map":
            # periodicidad de 90°: el mapa cubre 45°–135°
            php = (ph - 45.0) % 90.0 + 45.0
            from scipy.ndimage import map_coordinates
            ir = np.interp(r, self.r, np.arange(len(self.r)))
            ip = np.interp(php, self.ph, np.arange(len(self.ph)))
            Bz = map_coordinates(self.Bz, [ir, ip], order=1, mode="nearest")
            dr = self.r[1] - self.r[0]; dp = self.ph[1] - self.ph[0]
            dBdr = (map_coordinates(self.Bz, [ir + 1, ip], order=1, mode="nearest") -
                    map_coordinates(self.Bz, [ir - 1, ip], order=1, mode="nearest")) / (2 * dr)
            dBdp = (map_coordinates(self.Bz, [ir, ip + 1], order=1, mode="nearest") -
                    map_coordinates(self.Bz, [ir, ip - 1], order=1, mode="nearest")) / (2 * math.radians(dp))
        else:
            Bz = np.full_like(r, self.B0); dBdr = np.zeros_like(r); dBdp = np.zeros_like(r)
        Br = z * dBdr                                    # z en mm, dBdr en T/mm → T
        Bp = z * dBdp / r
        c, s = x / r, y / r
        return np.c_[Br * c - Bp * s, Br * s + Bp * c, Bz]


# ----------------------------------------------------------------------------- tracking
def track(Efun, Bfun, x0, pr, phi0_deg, t_max_turns=8, dt_ns=0.02, T0_eV=10.0, R_stop=130.0, z_stop=None):
    """Boris relativista. x0: (N,3) mm; dirección inicial = e_n (hacia el puller). Devuelve trayectorias."""
    q = -E0
    MM = pr.get("_mass", M_HMINUS_MEV)
    m = MM * 1e6 * E0 / C ** 2
    w = 2 * math.pi * pr["f_rf"]
    N = len(x0)
    ph0 = np.radians(np.broadcast_to(np.asarray(phi0_deg, float), (N,))).copy()
    x = np.array(x0, float) * 1e-3
    en = np.r_[np.array(pr["_e_n"]), 0.0]
    g0 = 1 + T0_eV / (MM * 1e6)
    p = np.tile(en * m * C * math.sqrt(g0 ** 2 - 1), (N, 1))
    dt = dt_ns * 1e-9
    nsteps = int(t_max_turns * pr["h"] / pr["f_rf"] / dt)
    alive = np.ones(N, bool)
    out = [[] for _ in range(N)]
    reasons = {}
    zs0 = z_stop or (pr["ap_half"] - 0.2)
    rcj = pr.get("r_caja", 0.0)                                  # v0.29: fuera de las puntas la abertura es la de la caja
    zs1 = min(pr.get("ap_caja", pr["ap_half"]), pr["hill_gap_half"]) - 0.2 if rcj else zs0
    t = 0.0
    for k in range(nsteps):
        idx = np.flatnonzero(alive)
        if not len(idx):
            break
        X = x[idx] * 1e3
        Eu = Efun(X)                                           # V/mm por volt de dee
        cosv = np.cos(w * t + ph0[idx])
        E = Eu * (pr["V_dee"] * cosv * 1e3)[:, None]           # V/m
        B = Bfun(X)
        pm = p[idx] + q * E * dt / 2
        gam = np.sqrt(1 + (np.einsum("ij,ij->i", pm, pm)) / (m * C) ** 2)
        tv = (q * dt / (2 * m * gam))[:, None] * B
        pp = pm + np.cross(pm, tv)
        sv = 2 * tv / (1 + np.einsum("ij,ij->i", tv, tv))[:, None]
        pm = pm + np.cross(pp, sv)
        p[idx] = pm + q * E * dt / 2
        gam = np.sqrt(1 + (np.einsum("ij,ij->i", p[idx], p[idx])) / (m * C) ** 2)
        x[idx] += p[idx] / (m * gam)[:, None] * dt
        t += dt
        Xn = x[idx] * 1e3
        r = np.hypot(Xn[:, 0], Xn[:, 1])
        metal = Efun.inside_metal(Xn)
        zs = np.where(r < rcj, zs0, zs1) if rcj else zs0
        bad = (r > R_stop) | (np.abs(Xn[:, 2]) > zs) | ~np.isfinite(Xn).all(axis=1) | metal
        if bad.any():
            for j in np.flatnonzero(bad):
                why = "r" if r[j] > R_stop else "z" if abs(Xn[j, 2]) > (zs[j] if np.ndim(zs) else zs) else "metal" if metal[j] else "nan"
                reasons.setdefault(int(idx[j]), why)
        if k % pr.get("_cada", 5) == 0 or bad.any():
            T = (gam - 1) * MM * 1e3                           # keV
            for j, i in enumerate(idx):
                out[i].append((t * 1e9, *Xn[j], T[j]))
        alive[idx[bad]] = False
    track.last_reasons = reasons
    return [np.array(o) for o in out]


class FEField:
    """E unitario de la malla FE (evaluación directa en el elemento) + detección de metal/fuera."""

    def __init__(self, mesh, gfu, mitad=False):
        from ngsolve import grad
        self.mesh = mesh
        self.E = -grad(gfu)
        self.mitad = mitad
        self._last_ok = None

    def _z(self, X):
        return np.maximum(np.abs(X[:, 2]), 1e-5) if self.mitad else X[:, 2]

    def __call__(self, X):
        mp = self.mesh(X[:, 0], X[:, 1], self._z(X))
        ok = mp["nr"] >= 0
        out = np.zeros((len(X), 3))
        if ok.any():
            out[ok] = np.asarray(self.E(mp[ok])).reshape(-1, 3)
        if self.mitad:
            out[:, 2] *= np.where(X[:, 2] < 0, -1.0, 1.0)          # espejo: E_z impar en z
        self._last_ok = ok
        return out

    def inside_metal(self, X):
        """Punto fuera de la malla de vacío → dentro de un electrodo o fuera del dominio (partícula perdida)."""
        mp = self.mesh(X[:, 0], X[:, 1], self._z(X))
        return mp["nr"] < 0


def solve(pr, out):
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, Preconditioner, TaskManager, \
        SetNumThreads, solvers
    t0 = time.time()
    air, info = build(pr)
    geo = OCCGeometry(air)
    ngm = geo.GenerateMesh(mp=mesh_params(pr, info))
    mesh = Mesh(ngm)
    print(f"[malla] {mesh.ne} tetraedros ({time.time() - t0:.0f} s); caras: {set(mesh.GetBoundaries())}", flush=True)
    SetNumThreads(2)
    fes = H1(mesh, order=2, dirichlet="tierra|dee|exterior")
    u, v = fes.TnT()
    gfu = GridFunction(fes)
    gfu.Set(1.0, definedon=mesh.Boundaries("dee"))
    a = BilinearForm(grad(u) * grad(v) * dx)
    if pr.get("_menisco") is not None:                          # v0.30: plasma detrás del menisco (φ = 0, penalización)
        a += 2e3 * chi_menisco(pr, info) * u * v * dx
    # bddc factoriza un problema grueso grande: con > 300 k tetraedros no cabe en 7 GB; h1amg es mucho más ligero
    import os
    pre = Preconditioner(a, os.environ.get("CR_PRECOND", "h1amg" if mesh.ne > 250000 else "bddc"))
    f = LinearForm(fes)
    with TaskManager():
        a.Assemble(); f.Assemble()
        r = f.vec.CreateVector(); r.data = f.vec - a.mat * gfu.vec
        inv = solvers.CGSolver(a.mat, pre.mat, maxiter=2000, tol=1e-10, printrates=False)
        du = gfu.vec.CreateVector(); du.data = inv * r
        # imponer solo en libres
        fd = np.array(fes.FreeDofs(), dtype=bool)
        gfu.vec.FV().NumPy()[fd] += du.FV().NumPy()[fd]
    print(f"[fem] {fes.ndof} gdl, {inv.iterations} it CG, {time.time() - t0:.0f} s", flush=True)
    return mesh, gfu, info, dict(ne=mesh.ne, ndof=fes.ndof, t_s=time.time() - t0)


# ----------------------------------------------------------------------------- análisis y figuras
def orbit_center(tr, turns_skip=2, pr=None):
    """Centro de giro medio a partir de las posiciones de la última vuelta (ajuste de círculo)."""
    if len(tr) < 20:
        return None
    T = tr[:, 0]; X = tr[:, 1]; Y = tr[:, 2]
    per = 1e9 * pr["h"] / pr["f_rf"]
    sel = T > T[-1] - per
    if sel.sum() < 10:
        return None
    A = np.c_[2 * X[sel], 2 * Y[sel], np.ones(sel.sum())]
    b = X[sel] ** 2 + Y[sel] ** 2
    cx, cy, c0 = np.linalg.lstsq(A, b, rcond=None)[0]
    return float(cx), float(cy), float(math.sqrt(max(c0 + cx * cx + cy * cy, 0)))


def _rf_markers(tr, pr, phi0_deg):
    """Posiciones donde V_dee = 0 (rojo) y |V_dee| máx (verde) a lo largo de una trayectoria."""
    if len(tr) < 3:
        return np.zeros((0, 2)), np.zeros((0, 2))
    w = 2 * math.pi * pr["f_rf"]
    ph = w * tr[:, 0] * 1e-9 + math.radians(phi0_deg)          # fase RF a lo largo de la órbita
    out = []
    for off in (math.pi / 2, 0.0):
        k0 = math.ceil((ph[0] - off) / math.pi); k1 = math.floor((ph[-1] - off) / math.pi)
        tgt = off + np.arange(k0, k1 + 1) * math.pi
        out.append(np.c_[np.interp(tgt, ph, tr[:, 1]), np.interp(tgt, ph, tr[:, 2])])
    return out[0], out[1]


def species_pr(pr, F):
    q = dict(pr)
    q.update(h=F["h"], _e_n=F["e_n"], _mass=ESPECIES[F["especie"]], _especie=F["especie"])
    return q


def run_species(Ef, pr, F, P, fases, vueltas, n_haz, R_stop, seed=3, beam=True):
    prs = species_pr(pr, F)
    Bf = MagField(prs, P)
    S0 = np.array(F["slit_center"]); en = np.array(F["e_n"])
    phases = np.arange(fases[0], fases[1] + 1e-9, fases[2])
    men = pr.get("_menisco")
    Sc = np.array(F["S"])
    if men is not None:                                       # v0.30: arranque sobre el menisco
        sm0 = men[0] + 0.05
        x0 = np.repeat([[Sc[0] + sm0 * en[0], Sc[1] + sm0 * en[1], 0.0]], len(phases), axis=0)
    else:
        x0 = np.repeat([[S0[0] + 0.2 * en[0], S0[1] + 0.2 * en[1], 0.0]], len(phases), axis=0)
    trs = track(Ef, Bf, x0, prs, phases, t_max_turns=vueltas, R_stop=R_stop)
    why = dict(track.last_reasons)
    t_full = 0.97 * vueltas * prs["h"] / prs["f_rf"] * 1e9
    res = []
    for i, (ph, tr) in enumerate(zip(phases, trs)):
        r_end = float(np.hypot(tr[-1, 1], tr[-1, 2])) if len(tr) else 0.0
        oc = orbit_center(tr, pr=prs)
        res.append(dict(fase=float(ph), T_final_keV=float(tr[-1, 4]) if len(tr) else 0.0, r_final_mm=r_end,
                        t_final_ns=float(tr[-1, 0]) if len(tr) else 0.0,
                        sobrevive=bool(len(tr) and (tr[-1, 0] > t_full or r_end > 0.97 * R_stop)),
                        centro_mm=oc, motivo=why.get(i, "vive")))
    ok = [r for r in res if r["sobrevive"]]
    dE_th = 4 * prs["V_dee"] * math.sin(math.radians(prs["h"] * prs["dee_angle"] / 2)) / 1e3
    best = None
    if ok:   # fase «de diseño»: la aceptada con el centro de órbita más cercano al centro de la máquina
        best = min(ok, key=lambda r: (math.hypot(r["centro_mm"][0], r["centro_mm"][1]) if r["centro_mm"] else 1e9))
    B = []
    if beam and ok and n_haz:
        rng = np.random.default_rng(seed)
        lo, hi = min(r["fase"] for r in ok), max(r["fase"] for r in ok)
        tv = np.array([-en[1], en[0]])
        PH = rng.uniform(lo - 10, hi + 10, n_haz) if pr.get("_PH") is None else np.asarray(pr["_PH"], float)
        nb = len(PH)
        U = rng.uniform(-0.45, 0.45, nb) * pr["slit_w"]; ZZ = rng.uniform(-0.4, 0.4, nb) * pr["slit_h"]
        if men is not None:
            sz = men[0] + men[1] * ZZ ** 2 + men[2] * ZZ ** 4 + men[3] * ZZ ** 6 + 0.05
            if len(men) > 4:
                sz = sz + men[4] * (2 * U / pr["slit_w"]) ** 2
            P0 = np.c_[Sc[0] + sz * en[0] + U * tv[0], Sc[1] + sz * en[1] + U * tv[1], ZZ]
        else:
            P0 = np.c_[S0[0] + 0.2 * en[0] + U * tv[0], S0[1] + 0.2 * en[1] + U * tv[1], ZZ]
        TR = track(Ef, Bf, P0, prs, PH, t_max_turns=vueltas, R_stop=R_stop)
        wb = dict(track.last_reasons)
        B = [dict(fase=float(ph), tr=tr, motivo=wb.get(k, "vive")) for k, (ph, tr) in enumerate(zip(PH, TR))]
    mk = None; long = None
    if best is not None:
        tb = trs[[r["fase"] for r in res].index(best["fase"])]
        red, green = _rf_markers(tb, prs, best["fase"])
        mk = dict(V0=red.tolist(), Vmax=green.tolist())
        if len(red):                                          # CAS 2017 Fig. 37: puntos V = 0 sobre los ejes de dees
            # V = 0 cada 180°/h de azimut: ejes de referencia = eje del dee + k·180°/h (h = 2: dees y dummy dees;
            # h = 4: además los centros de colina)
            ejes = np.radians(pr["dee_az"][0] + np.arange(2 * prs["h"]) * 180.0 / prs["h"])
            ang = np.arctan2(red[:, 1], red[:, 0])
            dev = np.degrees(np.min(np.abs((ang[:, None] - ejes[None, :] + np.pi) % (2 * np.pi) - np.pi), axis=1))
            long = dict(desviacion_media_deg=float(dev.mean()), desviacion_max_deg=float(dev.max()), n_puntos=int(len(dev)))
    vent_cas = None
    if ok:                                                    # guía CAS: fase de salida entre −40° y −10°
        acc = np.array([r["fase"] for r in ok])
        vent_cas = float(np.mean((acc >= -40) & (acc <= -10)))
    trans = float(np.mean([b["motivo"] in ("vive", "r") for b in B])) if B else 0.0   # "r": llegó al radio de salida
    return dict(especie=F["especie"], h=F["h"], campo_B=Bf.src, barrido=res,
                ventana_fase_aceptada=[min(r["fase"] for r in ok), max(r["fase"] for r in ok)] if ok else None,
                ancho_aceptancia_deg=(max(r["fase"] for r in ok) - min(r["fase"] for r in ok) + fases[2]) if ok else 0.0,
                dE_vuelta_gap_fino_keV=dE_th, fase_diseno=best, marcadores_RF=mk, ajuste_longitudinal=long,
                fraccion_ventana_en_guia_CAS=vent_cas,
                perdidas_haz={m_: sum(1 for b in B if b["motivo"] == m_) for m_ in ("vive", "z", "metal", "r", "nan")},
                transmision_haz=trans), trs, B, Bf


def optimize_sources(pr, P, radios_rel, fases=(-50.0, 10.0, 10.0), vueltas=3.0):
    """Optimización paramétrica: mueve las fuentes en radio (escala relativa) y elige, para cada especie,
    el radio cuyo centro de órbita (vuelta 3) queda más cerca del centro de la máquina."""
    hist = []
    base = [F["src_r"] for F in pr["fuentes"]]
    coarse = dict(h_gap=0.35, h_nose=1.6, h_mid=8.0, h_far=30.0, R_dom=min(pr["R_dom"], 220.0), R_fino=40.0)
    for k in radios_rel:
        q = dict(pr, **coarse)
        q["fuentes"] = [dict(F, src_r=b * k) for F, b in zip(pr["fuentes"], base)]
        mesh, gfu, info, st = solve(q, None)
        Ef = FEField(mesh, gfu, q.get("mitad", False))
        row = dict(escala=k, malla=st["ne"])
        for F in info["fuentes"]:
            r, *_ = run_species(Ef, q, F, P, fases, vueltas, 0, R_stop=0.9 * q["R_dom"], beam=False)
            ok = [x for x in r["barrido"] if x["sobrevive"] and x["centro_mm"]]
            dc = min((math.hypot(x["centro_mm"][0], x["centro_mm"][1]) for x in ok), default=None)
            row[F["especie"]] = dict(src_r=[G for G in q["fuentes"] if G["especie"] == F["especie"]][0]["src_r"],
                                     fases_ok=len(ok), desc_centro_mm=dc)
        print("[opt]", json.dumps(row), flush=True)
        hist.append(row)
    best = {}
    for F in pr["fuentes"]:
        cand = [h[F["especie"]] for h in hist if h[F["especie"]]["desc_centro_mm"] is not None]
        if cand:
            best[F["especie"]] = min(cand, key=lambda c: (c["desc_centro_mm"] - 0.5 * c["fases_ok"]))["src_r"]
    return best, hist


def export_electrodes(out, pr):
    """STEP de electrodos (para Opera-3D / Blender), recortados a la zona central para el render."""
    from netgen.occ import Cylinder, Pnt, Z
    d = out / "electrodos"; d.mkdir(exist_ok=True)
    clip = Cylinder(Pnt(0, 0, -100), Z, r=min(pr["R_dom"], 160.0), h=200)
    files = {}
    for name, shapes in build.groups.items():
        for i, sh in enumerate(shapes):
            f = d / f"{name}_{i}.step"
            try:
                (sh * clip).WriteStep(str(f)); files.setdefault(name, []).append(f.name)
            except Exception as e:                                # noqa: BLE001
                print("[aviso] STEP", name, i, e)
    return files


class GridE:
    """E unitario (por voltio de dee) leído de una tabla exportada de Opera-3d (u otro código):
    columnas X Y Z EX EY EZ en mm y V/mm, grilla regular. Permite rehacer el tracking con el campo de Opera."""

    def __init__(self, path, V_ref=1.0):
        from scipy.interpolate import RegularGridInterpolator
        d = np.loadtxt(path, comments=("#", "!"), delimiter=None if not str(path).endswith(".csv") else ",", skiprows=0)
        xs, ys, zs = (np.unique(d[:, k]) for k in range(3))
        sh = (len(xs), len(ys), len(zs))
        idx = np.lexsort((d[:, 2], d[:, 1], d[:, 0]))
        E = d[idx, 3:6].reshape(*sh, 3) / V_ref
        self.f = [RegularGridInterpolator((xs, ys, zs), E[..., k], bounds_error=False, fill_value=np.nan) for k in range(3)]

    def __call__(self, X):
        out = np.stack([f(X) for f in self.f], axis=1)
        self._nan = np.isnan(out).any(axis=1)
        return np.nan_to_num(out)

    def inside_metal(self, X):
        return np.isnan(np.stack([f(X) for f in self.f[:1]], axis=1)).any(axis=1)


def export_E_table(mesh, gfu, path, L=60.0, step=1.0, zs=(-6, -4.5, -3, -1.5, 0, 1.5, 3, 4.5, 6)):
    """Tabla X Y Z EX EY EZ (mm, V/mm por voltio de dee) del FE — mismo formato que se pide a Opera."""
    from ngsolve import grad
    E = -grad(gfu)
    xs = np.arange(-L, L + 1e-9, step)
    X, Y, Zz = np.meshgrid(xs, xs, np.asarray(zs, float), indexing="ij")
    pts = np.c_[X.ravel(), Y.ravel(), Zz.ravel()]
    za = np.maximum(np.abs(pts[:, 2]), 1e-5)                        # válido para malla completa o de mitad
    mp = mesh(pts[:, 0], pts[:, 1], za); ok = mp["nr"] >= 0
    val = np.full((len(pts), 3), np.nan)
    val[ok] = np.asarray(E(mp[ok])).reshape(-1, 3)
    val[:, 2] *= np.where(pts[:, 2] < 0, -1.0, 1.0)
    np.savetxt(path, np.c_[pts, val], fmt="%.4f %.4f %.4f %.6e %.6e %.6e",
               header="X[mm] Y[mm] Z[mm] EX EY EZ [V/mm por 1 V de dee]; NaN = dentro de un electrodo")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="kiube", choices=list(PRESETS))
    ap.add_argument("--set", nargs="*", default=[], help="clave=valor para modificar el preset (p. ej. gap_sp=3.0)")
    ap.add_argument("--fases", type=float, nargs=3, default=[-70.0, 30.0, 5.0], help="φ0 inicio fin paso [°RF]")
    ap.add_argument("--vueltas", type=float, default=6.0)
    ap.add_argument("--n_haz", type=int, default=120, help="partículas del haz por especie")
    ap.add_argument("--optimizar", type=float, nargs="*", default=None,
                    help="escalas relativas del radio de las fuentes a probar (p. ej. 0.8 0.9 1 1.1 1.2)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--campoE", default=None, help="tabla X Y Z EX EY EZ exportada de Opera (V/mm por voltio de dee)")
    ap.add_argument("--solo_campo", action="store_true", help="solo resolver y exportar la tabla de E (sin tracking)")
    a = ap.parse_args()
    pr = dict(PRESETS[a.preset]); pr["_nombre"] = a.preset
    for kv in a.set:
        k, v = kv.split("="); pr[k] = type(pr[k])(v) if not isinstance(pr[k], (list, dict)) else json.loads(v)
    out = pathlib.Path(a.out or ROOT / f"out/central_region/{a.preset}"); out.mkdir(parents=True, exist_ok=True)
    from kiube import load_params
    P, _ = load_params()
    opt = None
    if a.optimizar:
        best, hist = optimize_sources(pr, P, a.optimizar)
        pr["fuentes"] = [dict(F, src_r=best.get(F["especie"], F["src_r"])) for F in pr["fuentes"]]
        opt = dict(probado=hist, elegido=best)
        json.dump(opt, open(out / "optimizacion_fuentes.json", "w"), indent=1, default=float)
    mesh, gfu, info, st = solve(pr, out)
    export_E_table(mesh, gfu, out / "campo_E_referencia_FE.txt", L=min(60.0, pr["R_dom"]))
    if a.solo_campo:
        return
    Ef = FEField(mesh, gfu, pr.get("mitad", False)) if not a.campoE else GridE(a.campoE)
    R_stop = 0.9 * min(pr["R_dom"], 220.0)
    t0 = time.time()
    especies, tr_all, beams, Bfs = [], {}, {}, {}
    for F in info["fuentes"]:
        r, trs, B, Bf = run_species(Ef, pr, F, P, a.fases, a.vueltas, a.n_haz, R_stop)
        print(f"[{F['especie']}] {Bf.src}", flush=True)
        print(f"[{F['especie']}] ventana {r['ventana_fase_aceptada']}  diseño {r['fase_diseno'] and r['fase_diseno']['fase']}  "
              f"transmisión haz {r['transmision_haz']:.2f}", flush=True)
        especies.append(r); tr_all[F["especie"]] = trs; beams[F["especie"]] = B; Bfs[F["especie"]] = Bf
    summary = dict(preset=a.preset, parametros={k: v for k, v in pr.items() if not k.startswith("_")},
                   malla=st, fuentes=info["fuentes"], especies=especies, optimizacion=opt,
                   tiempo_tracking_s=time.time() - t0)
    json.dump(summary, open(out / "resultado.json", "w"), indent=1, default=float, ensure_ascii=False)
    npz = {}
    for r in especies:
        e = r["especie"].replace("-", "m")
        for x, t in zip(r["barrido"], tr_all[r["especie"]]):
            npz[f"{e}_fase_{x['fase']:+.1f}"] = t
        for k, b in enumerate(beams[r["especie"]]):
            npz[f"{e}_haz_{k:03d}"] = b["tr"]
    np.savez_compressed(out / "trayectorias.npz", **npz)
    try:
        summary["electrodos_step"] = export_electrodes(out, pr)
        json.dump(summary, open(out / "resultado.json", "w"), indent=1, default=float, ensure_ascii=False)
    except Exception as e:                                         # noqa: BLE001
        print("[aviso] exportación de electrodos:", e)
    figures_multi(pr, info, mesh, gfu, especies, tr_all, beams, out, Bfs)
    print(json.dumps([{k: r[k] for k in ("especie", "ventana_fase_aceptada", "dE_vuelta_gap_fino_keV", "transmision_haz")}
                      for r in especies], ensure_ascii=False))


COL = {"H-": "#1b7f3b", "D-": "#2a5db0"}


def _wedge(ax, R0, R1, a0, a1, **kw):
    from matplotlib.patches import Wedge
    ax.add_patch(Wedge((0, 0), R1, a0, a1, width=R1 - R0, **kw))


def _electrode_plan(ax, pr, info, L):
    from matplotlib.patches import Polygon, Circle
    for k in range(4):                                             # colinas (proyección)
        hc = 45.0 + 90 * k
        _wedge(ax, pr["plug_r"], L * 1.5, hc - pr["hill_angle"] / 2, hc + pr["hill_angle"] / 2,
               facecolor="#e7e9ee", edgecolor="#b5bac4", lw=0.6, zorder=0)
    ax.add_patch(Circle((0, 0), pr["plug_r"], facecolor="#eef0f3", edgecolor="#b5bac4", lw=0.6, zorder=0))
    for az in pr["dee_az"]:                                       # dees
        _wedge(ax, pr["nose_r"], L * 1.5, az - pr["dee_angle"] / 2, az + pr["dee_angle"] / 2,
               facecolor="#f3d6bd", edgecolor="#c0703a", lw=1.0, zorder=1, alpha=0.9)
        if pr.get("dummy"):
            _wedge(ax, pr["dummy_rmin"], L * 1.5, az + 90 - pr["dummy_angle"] / 2, az + 90 + pr["dummy_angle"] / 2,
                   facecolor="#dde3ec", edgecolor="#7d8899", lw=0.8, zorder=1, hatch="///", alpha=0.8)
    for F in info["fuentes"]:                                     # chimeneas, rendijas y pullers
        S, en = np.array(F["S"]), np.array(F["e_n"]); tv = np.array([-en[1], en[0]])
        ax.add_patch(Circle(S, pr["chimney_r"], facecolor="#5b6170", edgecolor="k", lw=0.5, zorder=4))
        c0 = S + (pr["chimney_r"] + pr["gap_sp"]) * en
        pts = [c0 + u * tv + v * en for u, v in ((-pr["puller_w"] / 2, 0), (pr["puller_w"] / 2, 0),
                                                 (pr["puller_w"] / 2, pr["puller_len"] + 9), (-pr["puller_w"] / 2, pr["puller_len"] + 9))]
        ax.add_patch(Polygon(pts, closed=True, facecolor="#e2a46f", edgecolor="#9c4f17", lw=0.6, zorder=3))
        ax.annotate(f"fuente {F['especie']}", S, xytext=(S[0] - 18 * en[0] + 10 * tv[0], S[1] - 18 * en[1] + 10 * tv[1]),
                    fontsize=8, color=COL.get(F["especie"], "k"), arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))


def figures_multi(pr, info, mesh, gfu, especies, tr_all, beams, out, Bfs):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    # ---- 1. planta tipo CAS Fig. 37: órbitas por especie, puntos V = 0 (rojo) y V máx (verde)
    rmax = max([np.hypot(t[:, 1], t[:, 2]).max() for e in especies for t in tr_all[e["especie"]] if len(t)] + [60])
    L = min(rmax * 1.08, pr["R_dom"])
    fig, ax = plt.subplots(figsize=(10, 10))
    _electrode_plan(ax, pr, info, L)
    for e in especies:
        c = COL.get(e["especie"], "k")
        for b in beams[e["especie"]][:60]:
            t = b["tr"]
            if len(t) and b["motivo"] == "vive":
                ax.plot(t[:, 1], t[:, 2], lw=0.3, color=c, alpha=0.25, zorder=5)
        d = e["fase_diseno"]
        if d is None:
            continue
        t = tr_all[e["especie"]][[x["fase"] for x in e["barrido"]].index(d["fase"])]
        ax.plot(t[:, 1], t[:, 2], lw=1.6, color=c, zorder=6,
                label=f"{e['especie']} (h = {e['h']}), φ₀ = {d['fase']:+.0f}°, aceptancia {e['ancho_aceptancia_deg']:.0f}° RF")
        mk = e["marcadores_RF"]
        if mk:
            r_, g_ = np.array(mk["V0"]), np.array(mk["Vmax"])
            ms_ = 6.0 if e["especie"] == "H-" else 4.2
            if len(r_):
                ax.plot(r_[:, 0], r_[:, 1], "o", ms=ms_, color="#d7191c", mec=c, mew=1.3, zorder=7)
            if len(g_):
                ax.plot(g_[:, 0], g_[:, 1], "o", ms=ms_, color="#4dd06a", mec=c, mew=1.3, zorder=7)
    ax.plot([], [], "o", color="#d7191c", label="V_dee = 0 (ideal: sobre los ejes de dees y dummy dees)")
    ax.plot([], [], "o", color="#4dd06a", label="V_dee máximo (borde = especie)")
    for az in pr["dee_az"]:
        ax.plot([0, L * math.cos(math.radians(az))], [0, L * math.sin(math.radians(az))], ls=":", lw=0.8, color="#c0703a")
    ax.set_aspect("equal"); ax.set_xlim(-L, L); ax.set_ylim(-L, L)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax.set_title(f"Región central {pr['_nombre']}: órbitas calculadas (E de FE 3D, Lorentz relativista, RF)\n"
                 f"f_RF {pr['f_rf'] / 1e6:.1f} MHz, V_dee {pr['V_dee'] / 1e3:.0f} kV, dees {pr['dee_angle']:.0f}°"
                 + (" + dummy dees" if pr.get('dummy') else ""), fontsize=10)
    fig.tight_layout(); fig.savefig(out / "orbitas_planta.png", dpi=160); plt.close(fig)
    # ---- 2. zoom fuente–puller con potencial en z = 0
    n = len(info["fuentes"])
    fig, axs = plt.subplots(1, n, figsize=(6.6 * n, 6.0), squeeze=False)
    for ax, F in zip(axs[0], info["fuentes"]):
        S = np.array(F["S"]); zl = 9.0
        xs = np.linspace(S[0] - zl, S[0] + zl, 320); ys = np.linspace(S[1] - zl, S[1] + zl, 320)
        X, Y = np.meshgrid(xs, ys)
        mp = mesh(X.ravel(), Y.ravel(), np.full(X.size, 1e-4)); ok = mp["nr"] >= 0
        phi = np.full(X.size, np.nan); phi[ok] = np.asarray(gfu(mp[ok])).ravel(); phi = phi.reshape(X.shape)
        ax.contourf(X, Y, np.where(np.isnan(phi), 1, np.nan), levels=[0.5, 1.5], colors=["#9aa0aa"])
        cs = ax.contourf(X, Y, phi, levels=np.linspace(0, 1, 21), cmap="coolwarm", alpha=0.85)
        for b in beams[F["especie"]][:50]:
            t = b["tr"]
            if len(t):
                ax.plot(t[:, 1], t[:, 2], lw=0.6, color="k", alpha=0.5)
        ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[0], ys[-1]); ax.set_aspect("equal")
        ax.set_title(f"Gap fuente–puller {F['especie']} (malla {pr['h_gap']} mm)", fontsize=9)
    fig.colorbar(cs, ax=axs[0].tolist(), label="φ / V_dee", shrink=0.85)
    fig.savefig(out / "fuente_puller_zoom.png", dpi=150); plt.close(fig)
    # ---- 3. aceptancia en fase y energía
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for e in especies:
        c = COL.get(e["especie"], "k")
        f_ = [x["fase"] for x in e["barrido"]]
        ax[0].plot(f_, [x["T_final_keV"] if x["sobrevive"] else np.nan for x in e["barrido"]], "o-", color=c,
                   label=f"{e['especie']} sobrevive")
        ax[0].plot(f_, [x["T_final_keV"] if not x["sobrevive"] else np.nan for x in e["barrido"]], "x", color=c, alpha=0.5)
        d = e["fase_diseno"]
        if d:
            t = tr_all[e["especie"]][f_.index(d["fase"])]
            Tturn = 1e9 * e["h"] / pr["f_rf"]
            ax[1].plot(t[:, 0] / Tturn, t[:, 4], color=c, lw=1.4, label=f"{e['especie']} φ₀ = {d['fase']:+.0f}°")
            nt = np.linspace(0, t[-1, 0] / Tturn, 50)
            ax[1].plot(nt, nt * e["dE_vuelta_gap_fino_keV"], ls="--", color=c, lw=0.8,
                       label=f"{e['especie']}: máximo de gap fino 4·V·sin(hα/2) = {e['dE_vuelta_gap_fino_keV']:.0f} keV/vuelta")
    ax[0].set_xlabel("fase RF de salida φ₀ [°]"); ax[0].set_ylabel("T al final [keV]"); ax[0].set_title("Aceptancia en fase (× = perdida)")
    ax[1].set_xlabel("vueltas"); ax[1].set_ylabel("T [keV]"); ax[1].set_title("Ganancia de energía")
    for x_ in ax:
        x_.grid(alpha=.3); x_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "aceptancia_energia.png", dpi=140); plt.close(fig)
    # ---- 4. movimiento vertical
    fig, ax = plt.subplots(figsize=(9, 3.8))
    for e in especies:
        for b in beams[e["especie"]][:60]:
            t = b["tr"]
            if len(t):
                ax.plot(np.hypot(t[:, 1], t[:, 2]), t[:, 3], lw=0.4, color=COL.get(e["especie"], "k"), alpha=0.6)
    ax.axhline(pr["ap_half"], color="#c0703a", lw=0.8, ls="--"); ax.axhline(-pr["ap_half"], color="#c0703a", lw=0.8, ls="--")
    ax.set_xlabel("r [mm]"); ax.set_ylabel("z [mm]"); ax.set_title("Movimiento vertical (líneas: abertura del dee)")
    ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(out / "vertical.png", dpi=140); plt.close(fig)
    # ---- 5. fuerza vertical eléctrica a lo largo de la órbita de diseño (equivalente a CAS 2017 Fig. 40)
    Ef = FEField(mesh, gfu, pr.get("mitad", False))
    fig, axs = plt.subplots(len(especies), 1, figsize=(11, 3.2 * len(especies)), squeeze=False)
    for ax, e in zip(axs[:, 0], especies):
        d = e["fase_diseno"]
        if d is None:
            continue
        t = tr_all[e["especie"]][[x["fase"] for x in e["barrido"]].index(d["fase"])]
        X = np.c_[t[:, 1], t[:, 2], np.full(len(t), 1.0)]                # sonda a z = +1 mm del plano medio
        Ez = Ef(X)[:, 2]
        V = np.cos(2 * math.pi * pr["f_rf"] * t[:, 0] * 1e-9 + math.radians(d["fase"]))
        Fz = -1.0 * Ez * V                                              # q < 0: F_z = q·E_z·V(t) (signo)
        Fz /= np.max(np.abs(Fz)) or 1.0
        Tt = 1e9 * e["h"] / pr["f_rf"]
        ax.plot(t[:, 0] / Tt, Fz, color=COL.get(e["especie"], "k"), lw=1)
        ax.axhline(0, color="0.5", lw=0.6)
        ph = np.degrees(np.arctan2(t[:, 2], t[:, 1])) % 360
        for az in pr["dee_az"]:                                       # cruces de los bordes de los dees = gaps
            for edge in (az - pr["dee_angle"] / 2, az + pr["dee_angle"] / 2):
                dd = ((ph - edge + 180) % 360) - 180
                k = np.flatnonzero((dd[:-1] < 0) & (dd[1:] >= 0) & (np.abs(dd[:-1]) < 20))
                for i in k:
                    ax.axvline(t[i, 0] / Tt, color="#2e8b57", lw=0.5, alpha=0.7)
        ax.set_xlim(0, min(5, t[-1, 0] / Tt)); ax.set_ylabel("F_z normalizada (z = +1 mm)")
        ax.set_title(f"{e['especie']}: fuerza vertical eléctrica (negativo = enfoca hacia el plano medio); verde = cruces de gap", fontsize=9)
        ax.grid(alpha=.3)
    axs[-1, 0].set_xlabel("vueltas")
    fig.tight_layout(); fig.savefig(out / "fuerza_vertical.png", dpi=140); plt.close(fig)


if __name__ == "__main__":
    main()
