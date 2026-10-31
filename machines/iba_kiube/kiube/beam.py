"""Óptica hard-edge 2D del KIUBE: órbita cerrada a 18 MeV y trayectorias de extracción.

v0.2: campos calibrados con el mapa medido (SUPAG09 Fig.1), ancho efectivo de colina
resuelto por isocronía exacta con f_RF = 40.646 MHz (TUP04 Fig.8), h = 2, y strippers
ajustados para que el haz extraído pase por el centro de los puertos medidos en
TUP04 Fig.2. Modelo simplificado: NO sustituye un mapa OPERA/AOC.
"""
import math
import numpy as np

C = 299.792458          # MeV/(T·m)
M_HMINUS = 939.294      # MeV (p + 2e)
CLIGHT = 299792458.0


def brho(T, m=M_HMINUS):
    return math.sqrt(T * T + 2 * T * m) / C


def wrap(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


class HardEdgeField:
    def __init__(self, P, hill_angle=None):
        p, b = P["pole"], P["beam"]
        self.n = p["n_sectors"]
        self.hw = math.radians(hill_angle if hill_angle else b.get("hill_eff_angle", p["hill_angle"])) / 2
        self.c0 = math.radians(p["hill_center0"])
        self.Bh, self.Bv = b["B_hill_T"], b["B_valley_T"]
        self.edge = None          # (d0, R) mm→m, arco del escalón 2
        self.r_valley = 9.0

    def B(self, x, y):
        th = math.atan2(y, x)
        per = 2 * math.pi / self.n
        k = round((th - self.c0) / per)
        d = wrap(th - (self.c0 + k * per))
        if abs(d) <= self.hw:
            if self.edge is None:
                return self.Bh
            d0, R = self.edge
            c = self.c0 + k * per
            return self.Bh if math.hypot(x - d0 * math.cos(c), y - d0 * math.sin(c)) <= R else 0.0
        return self.Bv if math.hypot(x, y) <= self.r_valley else 0.0


def track(F, x, y, ang, sign, br, ds=2e-4, stop=None, nmax=200000):
    pts = [(x, y, ang)]
    for _ in range(nmax):
        k = sign * F.B(x, y) / br
        dang = k * ds
        if abs(dang) > 1e-12:
            R = 1 / k
            x += R * (math.sin(ang + dang) - math.sin(ang))
            y -= R * (math.cos(ang + dang) - math.cos(ang))
        else:
            x += ds * math.cos(ang); y += ds * math.sin(ang)
        ang += dang
        pts.append((x, y, ang))
        if stop and stop(x, y, ang, pts):
            break
    return np.array(pts)


def closed_orbit(F, br):
    c, per = F.c0, 2 * math.pi / F.n

    def err(r0):
        x, y = r0 * math.cos(c), r0 * math.sin(c)

        def stop(x, y, a, pts):
            return wrap(math.atan2(y, x) - (c + per)) >= 0 and len(pts) > 10
        t = track(F, x, y, c + math.pi / 2, +1, br, stop=stop)
        xe, ye, ae = t[-1]
        return wrap(ae - (math.atan2(ye, xe) + math.pi / 2)), t

    lo, hi = 0.30, 0.65
    elo, _ = err(lo)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        em, t = err(mid)
        if np.sign(em) == np.sign(elo):
            lo, elo = mid, em
        else:
            hi = mid
    return mid, t


def rev_freq(orb_quarter, T, n=4, m=M_HMINUS):
    xy = orb_quarter[:, :2]
    L = np.sum(np.hypot(*np.diff(xy, axis=0).T)) * n
    g = 1 + T / m
    return math.sqrt(1 - 1 / g ** 2) * CLIGHT / L


def solve_isochronous_hill(P):
    """Ancho efectivo de colina tal que f_rev(18 MeV) = f_RF / h."""
    b = P["beam"]
    br = brho(b["T_MeV"])
    f_target = b["rf_MHz"] * 1e6 / b["rf_harmonic"]
    lo, hi = 35.0, 65.0
    for _ in range(22):
        m = 0.5 * (lo + hi)
        F = HardEdgeField(P, hill_angle=m)
        _, o = closed_orbit(F, br)
        if rev_freq(o, b["T_MeV"]) < f_target:
            lo = m
        else:
            hi = m
    return m


def de_voltage_MeV(b):
    return b["dee_voltage_kV"] / 1000.0


def cruce_cara_yugo(x, y, a, H):
    """Punto donde la recta (x, y) + s·(cos a, sin a), s > 0, sale del cuadrado |x|, |y| ≤ H (cara del yugo) [m]."""
    c, sn = math.cos(a), math.sin(a)
    ss = [s_ for s_ in ((H - x) / c if c > 1e-12 else None, (-H - x) / c if c < -1e-12 else None,
                        (H - y) / sn if sn > 1e-12 else None, (-H - y) / sn if sn < -1e-12 else None) if s_ is not None and s_ > 0]
    s_ = min(ss)
    return x + s_ * c, y + s_ * sn


def ventana_blanco(P, t):
    """Posición y dirección de la ventana del blanco a partir del final de la trayectoria t (fuera del campo: recta).
    v0.27: en la cara del yugo + blanco_off_cara (VE24 Fig. 3: los blancos se montan fuera de la cara, a ±0.87 m;
    foto VUB 2024: puertas blindadas en la cara). Si extraction.blanco_en_cara es falso, se usa window_radius."""
    ex, Y = P["extraction"], P["yoke"]
    x, y, a = t[-1, 0], t[-1, 1], t[-1, 2]
    if ex.get("blanco_en_cara", True):
        xc, yc = cruce_cara_yugo(x, y, a, Y["side"] / 2000)
        off = ex.get("blanco_off_cara", 10.0) / 1000
        return xc + off * math.cos(a), yc + off * math.sin(a), a
    return None


def extraction(P):
    b, po, ex = P["beam"], P["pole"], P["extraction"]
    if b.get("hill_eff_angle") in (None, "auto"):
        b["hill_eff_angle"] = solve_isochronous_hill(P)
    F = HardEdgeField(P)
    br = brho(b["T_MeV"])
    r0, orb1 = closed_orbit(F, br)
    f_rev = rev_freq(orb1, b["T_MeV"])
    rho_h = br / F.Bh
    L = po["step2_edge_L"] / 1000          # medido (TUP04 Fig.2)
    clr = L - r0
    d0 = r0 - rho_h
    F.edge = (d0, rho_h + clr)
    F.r_valley = L
    per = 2 * math.pi / F.n
    xy = orb1[:, :2]
    orbit = np.vstack([np.c_[xy[:, 0] * math.cos(k * per) - xy[:, 1] * math.sin(k * per),
                             xy[:, 0] * math.sin(k * per) + xy[:, 1] * math.cos(k * per)]
                       for k in range(F.n)])
    oang = np.arctan2(orbit[:, 1], orbit[:, 0])

    Rw = ex["window_radius"] / 1000
    Rp = ex["port_ref_radius"] / 1000

    def shoot(phs):
        i = int(np.argmin(np.abs((oang - phs + np.pi) % (2 * np.pi) - np.pi)))
        j = (i + 1) % len(orbit)
        x, y = orbit[i]
        a = math.atan2(orbit[j, 1] - y, orbit[j, 0] - x)
        st = {}

        def stop(x_, y_, a_, pts):
            if "edge" not in st and F.B(x_, y_) == 0.0 and math.hypot(x_, y_) > 0.3:
                st["edge"] = (x_, y_)
            if "win" not in st and math.hypot(x_, y_) >= Rw:
                st["win"] = (x_, y_, a_)
            return math.hypot(x_, y_) >= Rp
        t = track(F, x, y, a, -1, br, stop=stop)
        return (x, y), st, t

    strippers, targets, edges, paths, ports = [], [], [], [], []
    for k in range(F.n):
        hc = F.c0 + k * per
        for pk, poff in enumerate(ex["port_offsets"]):
            want = hc + math.radians(poff)
            lo, hi = hc - math.radians(25), hc + math.radians(22)

            def f(o):
                _, _, t = shoot(o)
                xe, ye = t[-1, :2]
                return wrap(math.atan2(ye, xe) - want)
            flo = f(lo)
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                fm = f(mid)
                if np.sign(fm) == np.sign(flo):
                    lo, flo = mid, fm
                else:
                    hi = mid
            (sx, sy), st, t = shoot(mid)
            wx, wy, wa = ventana_blanco(P, t) or st["win"]
            pxy = t[-1]
            strippers.append(dict(x=sx * 1e3, y=sy * 1e3, phi=math.degrees(mid),
                                  dphi_hill=math.degrees(wrap(mid - hc)), pole=k, sys=pk))
            targets.append(dict(x=wx * 1e3, y=wy * 1e3, dir=math.degrees(wa),
                                phi=math.degrees(math.atan2(wy, wx)), pole=k, sys=pk))
            ports.append(dict(x=pxy[0] * 1e3, y=pxy[1] * 1e3, dir=math.degrees(pxy[2]),
                              phi=math.degrees(math.atan2(pxy[1], pxy[0])), pole=k, sys=pk))
            e = st.get("edge", (st["win"][0], st["win"][1]) if "win" in st else (wx, wy))
            edges.append(dict(x=e[0] * 1e3, y=e[1] * 1e3, phi=math.degrees(math.atan2(e[1], e[0])),
                              dphi_hill=math.degrees(wrap(math.atan2(e[1], e[0]) - hc)), pole=k, sys=pk))
            paths.append(t[:, :2] * 1e3)
    # ganancia por vuelta y separación de vueltas (criterio de posición radial del stripper,
    # Kazarinov et al. TR24, Cyclotrons2019: R_f = R_ex − ΔR/2)
    de = P["dee"]
    dE = 2 * de["n"] * de_voltage_MeV(b) * math.sin(math.radians(b["rf_harmonic"] * de["angle"] / 2))
    if b.get("dE_turn_keV") in (None, "auto"):
        b["dE_turn_keV"] = dE * 1e3
    T = b["T_MeV"]
    pc2 = T * T + 2 * T * M_HMINUS
    dR = r0 * 1e3 * (T + M_HMINUS) * (b["dE_turn_keV"] / 1e3) / pc2
    res = dict(dE_turn_keV=b["dE_turn_keV"], dR_turn_mm=dR, n_turns=T * 1e3 / b["dE_turn_keV"],
                r_ext_hill=r0 * 1e3, r_ext_valley=float(np.hypot(*xy.T).min() * 1e3),
                f_rev_MHz=f_rev / 1e6, hill_eff_angle=b["hill_eff_angle"],
                brho=br, rho_hill=rho_h * 1e3, edge_d0=d0 * 1e3, edge_R=(rho_h + clr) * 1e3,
                edge_L=L * 1e3, Rh_over_Lh=(rho_h + clr) / L,
                orbit=orbit * 1e3, strippers=strippers, targets=targets, ports=ports,
                edge_crossings=edges, paths=paths)
    return apply_variable_energy(P, res)


def variable_energy(P, E, T_MeV, sel=None):
    """Energía variable (do Carmo et al., Instruments 8 (2024) 38 [VE24]).

    Un stripper de horquilla larga corta el haz en la órbita de energía T < 18 MeV; el protón debe
    salir por la MISMA ventana de la cámara (válvula y puerto fijos) y una brida en cuña alinea el
    blanco con la nueva dirección. Para cada puerto: azimut del stripper que lleva el haz al punto
    de la ventana de 18 MeV, radio de stripping, ángulo de cuña y cruce a r = 820 mm.
    sel: lista de índices de puerto (0..7) o None = todos.
    """
    b, ex = P["beam"], P["extraction"]
    F = HardEdgeField(P)
    F.edge = (E["edge_d0"] / 1e3, E["edge_R"] / 1e3)
    F.r_valley = E["edge_L"] / 1e3
    br = brho(T_MeV)
    r0, orb1 = closed_orbit(F, br)
    per = 2 * math.pi / F.n
    xy = orb1[:, :2]
    orbit = np.vstack([np.c_[xy[:, 0] * math.cos(k * per) - xy[:, 1] * math.sin(k * per),
                             xy[:, 0] * math.sin(k * per) + xy[:, 1] * math.cos(k * per)] for k in range(F.n)])
    oang = np.arctan2(orbit[:, 1], orbit[:, 0])
    Rw, Rp = ex["window_radius"] / 1e3, ex["port_ref_radius"] / 1e3

    def shoot(phs):
        i = int(np.argmin(np.abs((oang - phs + np.pi) % (2 * np.pi) - np.pi)))
        j = (i + 1) % len(orbit)
        x, y = orbit[i]
        a = math.atan2(orbit[j, 1] - y, orbit[j, 0] - x)
        st = {}

        def stop(x_, y_, a_, pts):
            if "win" not in st and math.hypot(x_, y_) >= Rw:
                st["win"] = (x_, y_, a_)
            return math.hypot(x_, y_) >= Rp + 0.05
        t = track(F, x, y, a, -1, br, stop=stop)
        return (x, y), st, t

    def cross(t, R):
        r = np.hypot(t[:, 0], t[:, 1]); i = int(np.argmax(r >= R))
        return math.degrees(math.atan2(t[i, 1], t[i, 0]))

    out = []
    for idx, (s18, w18) in enumerate(zip(E["strippers"], E["targets"])):
        if sel is not None and idx not in sel:
            continue
        hc = F.c0 + s18["pole"] * per
        want = math.atan2(w18["y"], w18["x"])
        lo, hi = hc - math.radians(70), hc + math.radians(25)

        def win_of(st, t):
            vb = ventana_blanco(P, t)
            return vb if vb else st.get("win")

        def f(o):
            _, st, t = shoot(o)
            w = win_of(st, t)
            if w is None:
                return -1.0
            return wrap(math.atan2(w[1], w[0]) - want)
        flo = f(lo)
        for _ in range(32):
            mid = 0.5 * (lo + hi); fm = f(mid)
            if np.sign(fm) == np.sign(flo):
                lo, flo = mid, fm
            else:
                hi = mid
        (sx, sy), st, t = shoot(mid)
        wx, wy, wa = win_of(st, t)
        out.append(dict(port=idx, pole=s18["pole"], sys=s18["sys"], T_MeV=T_MeV,
                        r_strip=math.hypot(sx, sy) * 1e3, phi_strip=math.degrees(mid),
                        dphi_hill_strip=math.degrees(wrap(mid - hc)),
                        win_miss_mm=math.hypot(wx - w18["x"] / 1e3, wy - w18["y"] / 1e3) * 1e3,
                        wedge_deg=math.degrees(wrap(wa - math.radians(w18["dir"]))),
                        phi_700=cross(t, 0.70) % 360, phi_820=cross(t, Rp) % 360,
                        dphi_hill_820=math.degrees(wrap(math.radians(cross(t, Rp)) - hc)),
                        path=t[:, :2] * 1e3))
    return dict(T_MeV=T_MeV, r_orbit_hill=r0 * 1e3, ports=out)


def apply_variable_energy(P, E):
    """Si extraction.ve_energy_MeV < 18: mueve los strippers de los puertos VE a la horquilla larga y
    gira el eje de sus blancos por el ángulo de la cuña (el punto de la ventana no cambia) [VE24]."""
    ex = P["extraction"]
    T = float(ex.get("ve_energy_MeV", P["beam"]["T_MeV"]))
    sel = list(ex.get("ve_ports", []))
    if T >= P["beam"]["T_MeV"] - 1e-6 or not sel:
        E["ve"] = None
        return E
    cunas = ex.get("ve_energias_cuna")
    if cunas and not any(abs(T - c) < 1e-6 for c in cunas):
        import warnings
        warnings.warn(f"ve_energy_MeV = {T}: IBA solo define cuñas para {cunas} MeV [VE24]; se calcula igual")
    V = variable_energy(P, E, T, sel=sel)
    for p in V["ports"]:
        i = p["port"]
        a = math.radians(p["phi_strip"])
        s = E["strippers"][i]
        s.update(x=p["r_strip"] * math.cos(a), y=p["r_strip"] * math.sin(a), phi=p["phi_strip"],
                 dphi_hill=p["dphi_hill_strip"], long_fork_mm=math.hypot(s["x"], s["y"]) - p["r_strip"])
        E["targets"][i]["dir"] += p["wedge_deg"]
        E["targets"][i]["wedge_deg"] = p["wedge_deg"]
        E["paths"][i] = p["path"]
    E["ve"] = dict(T_MeV=T, ports=[{k: v for k, v in p.items() if k != "path"} for p in V["ports"]])
    return E
