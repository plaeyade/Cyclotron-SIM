"""Órbitas en la región central del KIUBE con física real (Fase 68).

Estado: SIMULACIÓN TEÓRICA (memoria 3: sin haz real; validación teórica).
No sustituye OPERA/AOC; es tracking propio con anclas MEDIDAS y RF real.

Campo magnético (modelo sector hard-edge con perfil radial):
- Anclas medidas S09 (Van de Walle et al., ICAP2018 SUPAG09 Fig.1):
  colina 2,2–2,3 T → B_colina(r<300 mm)=2,25 T; sobre inserto (r≥480 mm)
  1,8–1,9 T → 1,85 T; interpolación lineal entre ambas (CALC declarado).
  Valle 0,40–0,50 T a r=350–500 mm → B_valle=0,45 T.
- Magnetostática propia: B ∝ 1/gap. Chequeo: B_colina·gap_colina/gap_valle
  = 2,15·24/132 = 0,39 T vs 0,45 T medido (desviación 13 %, borde/fringing;
  declarado). Se usan los valores MEDIDOS, no el escalado.
- Más allá de r_valle_ef (escalón 2, T04F2: L=512 mm) el campo decae a 0
  (borde de colina, arco descentrado — PAT).

RF real:
- f_RF = 40,646 MHz (T04F8), h = 2 (CALC: B0=1,335 T), V_dee = 40 kV (TUD03).
- 2 dees en valles 90°/270° (TUD03), abertura α_dee = 38° (EST declarado).
- Ganancia por cruce de gap: ΔE = V·sin(h·α/2)·cos(φ_RF(t_cruce)) — fase
  evaluada con el TIEMPO ABSOLUTO de vuelo (no sincronismo forzado).

Integrador: paso circular exacto en B local (mismo esquema que kiube v0.4
beam.py::track, verificado) con acumulación de tiempo ds/v y kicks de RF.

Isocronía: ancho efectivo de colina resuelto para que f_rev(18 MeV) =
f_RF/h (CALC, mismo criterio que kiube v0.4 → ≈46,7°).

Fuentes (49.11): TUP03, TUP04/T04F2/T04F8, TUD03, S09, PAT, CALC, EST.
PQ: no_declarada_nunca (14.9). Recomendación ≠ decisión operacional (18.10).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

C_MEV_TM = 299.792458          # MeV/(T·m)
M_HMINUS = 939.294             # MeV/c² (p + 2e)
CLIGHT = 299792458.0           # m/s

FUENTES_ORBITA = {
    "S09": "Van de Walle et al., ICAP2018 SUPAG09 Fig.1 — mapa de campo MEDIDO "
           "(colina 2,2–2,3 T; inserto 1,8–1,9 T; valle 0,40–0,50 T)",
    "T04F8": "Zaremba et al., Cyclotrons2016 TUP04 Fig.8 — f_RF = 40,646 MHz",
    "TUD03": "Nactergal et al., Cyclotrons2016 — V_dee = 40 kV; dees en valles 90°/270°",
    "TUP03": "Kleeven et al., Cyclotrons2016 — T_final = 18 MeV, ion H-",
    "T04F2": "TUP04 Fig.2 (planta a escala 937 px = 1740 mm) — L escalón 2 = 512 mm",
    "CALC": "Derivado por cálculo de datos publicados (isocronía, perfiles)",
    "EST": "Estimado — reemplazar cuando haya documento",
    "MAG-1GAP": "Magnetostática propia B∝1/gap: chequeo valle 0,39 T vs 0,45 T medido (13 %, declarado)",
    "KEY22": "Nuttens et al., Cyclotrons2022 WEAO03 Fig.3 — ⟨B⟩(r) y flutter "
             "F(r) del KIUBE MEDIDOS/calculados por IBA (digitalizado, "
             "±15 G, ±0,005 en F; F ambiguo en 0,285–0,345 m, interpolado)",
    "KEY22-GEOM": "KEY22 Fig.3 fijó por FE 3D (kiube v0.7): valle a 110 mm del "
                  "plano medio (perfil de flutter coincide) y tapón rebajado "
                  "17 mm (joroba central 1,04 vs 1,03 medida)",
}


def brho_mev_tm(T_mev: float, m: float = M_HMINUS) -> float:
    return math.sqrt(T_mev * T_mev + 2 * T_mev * m) / C_MEV_TM


def wrap_pi(a: float) -> float:
    return (a + math.pi) % (2 * math.pi) - math.pi


@dataclass
class CampoKIUBE:
    """Campo sector hard-edge isócrono anclado a medidas S09 (magnetostática propia).

    Modelo: para cada radio se impone isocronía exacta  B̄(r) = γ(r)·B0  con
    γ(r) = 1/√(1−(0,2998·B0·r/0,9393)²)  (órbita de equilibrio, CALC), y la
    razón medida ρ = B_valle/B_colina = 0,45/2,15 = 0,2093 (S09). Con el
    ancho efectivo de colina w (resuelto anclando B_colina(300 mm) = 2,25 T,
    S09) resultan B_colina(r) y B_valle(r) punto a punto. Así el campo es
    isócrono POR CONSTRUCCIÓN y honra las dos anclas medidas; no es un mapa
    medido completo (pendiente declarado).
    """
    n_sectores: int = 4                  # TUP04
    colina_centro0_deg: float = 45.0     # T04F2
    colina_ancho_ef_deg: float = 44.3    # CALC (ancla S09; re-soluble)
    b0_t: float = 1.3352                 # CALC: f_RF/h a γ=1 (B0 isócrono)
    rho_vh: float = 0.45 / 2.15          # S09 medido: B_valle/B_colina
    bh_ancla_t: float = 2.25             # S09 medido (2,2–2,3) en colina
    r_ancla_mm: float = 300.0            # CALC (radio del ancla)
    r_valle_ef_mm: float = 512.0         # T04F2 (escalón 2)
    r_camara_mm: float = 548.0           # T04F2
    shim_b0: float = 1.0                 # CALC: shim isócrono (calibrar_shim)
    shim_r_mm: tuple = ()                # CALC: shim radial (tabla)
    shim_f: tuple = ()
    # F69: perfil MEDIDO KEY22 (⟨B⟩(r) y flutter F(r)). Si se cargan, b_media
    # interpola la medida (×shim) en vez de γ·B0, y b_colina/b_valle salen del
    # flutter medido: F = ⟨B²⟩/⟨B⟩²−1 = a(1−a)(Bh−Bv)²/B̄² con a = w/90.
    key22_r_mm: tuple = ()               # KEY22 (medido, digitalizado)
    key22_bavg_t: tuple = ()             # KEY22
    key22_f: tuple = ()                  # KEY22 (huecos interpolados, declarado)
    gap_valle_mm: float = 132.0          # 110 con KEY22-GEOM (F69)

    def gamma_orbita(self, r_mm: float) -> float:
        a = 0.299792458 * self.b0_t * self.shim_b0 * (r_mm / 1e3) / 0.939294
        return 1.0 / math.sqrt(max(1 - a * a, 1e-12))

    def b_media(self, r_mm: float) -> float:
        if self.key22_bavg_t:
            b = float(np.interp(r_mm, self.key22_r_mm, self.key22_bavg_t))
            b *= self.shim_b0
        else:
            b = self.gamma_orbita(r_mm) * self.b0_t * self.shim_b0
        if self.shim_r_mm:
            b *= float(np.interp(r_mm, self.shim_r_mm, self.shim_f))
        return b

    def b_colina(self, r_mm: float) -> float:
        w = self.colina_ancho_ef_deg
        if self.key22_f:
            a = w / 90.0
            bm = self.b_media(r_mm)
            f = float(np.interp(r_mm, self.key22_r_mm, self.key22_f))
            d = bm * math.sqrt(max(f, 0.0) / (a * (1 - a)))
            return bm + (1 - a) * d
        return self.b_media(r_mm) * 90.0 / (w + self.rho_vh * (90.0 - w))

    def b_valle(self, r_mm: float) -> float:
        if self.key22_f:
            a = self.colina_ancho_ef_deg / 90.0
            bm = self.b_media(r_mm)
            f = float(np.interp(r_mm, self.key22_r_mm, self.key22_f))
            d = bm * math.sqrt(max(f, 0.0) / (a * (1 - a)))
            return bm - a * d
        return self.rho_vh * self.b_colina(r_mm)

    def resolver_ancho_colina(self) -> float:
        """w tal que B_colina(r_ancla) = bh_ancla (ancla S09). CALC."""
        w = (self.b_media(self.r_ancla_mm) * 90.0 / self.bh_ancla_t
             - self.rho_vh * 90.0) / (1 - self.rho_vh)
        self.colina_ancho_ef_deg = w
        return w

    def B(self, x_m: float, y_m: float) -> float:
        r_mm = math.hypot(x_m, y_m) * 1e3
        if r_mm > self.r_valle_ef_mm:
            return 0.0
        th = math.atan2(y_m, x_m)
        per = 2 * math.pi / self.n_sectores
        c0 = math.radians(self.colina_centro0_deg)
        k = round((th - c0) / per)
        d = wrap_pi(th - (c0 + k * per))
        if abs(d) <= math.radians(self.colina_ancho_ef_deg) / 2:
            return self.b_colina(r_mm)
        return self.b_valle(r_mm)


@dataclass
class RFKIUBE:
    f_mhz: float = 40.646                # T04F8
    armonico: int = 2                    # CALC
    v_dee_kv: float = 40.0               # TUD03
    n_dee: int = 2                       # TUD03 (valles 90/270)
    alfa_dee_deg: float = 38.0           # EST declarado
    valles_deg: tuple = (90.0, 270.0)    # TUD03

    @property
    def omega(self) -> float:
        return 2 * math.pi * self.f_mhz * 1e6

    @property
    def ganancia_gap_kev(self) -> float:
        """Referencia analítica: ganancia pico por cruce de borde de dee."""
        return (self.v_dee_kv
                * math.sin(math.radians(self.armonico * self.alfa_dee_deg / 2)))

    def v_dee_instantaneo_kv(self, t_s: float) -> float:
        """Potencial del dee (ambos en fase, h=2 — ver docstring del módulo)."""
        return self.v_dee_kv * math.cos(self.omega * t_s)

    def lineas_gap(self) -> list[tuple[float, int, int]]:
        """(azimut_rad, dee, dirección): dirección +1 = entrando al dee."""
        ls = []
        for i, v in enumerate(self.valles_deg):
            ls.append((math.radians(v - self.alfa_dee_deg / 2), i, +1))
            ls.append((math.radians(v + self.alfa_dee_deg / 2), i, -1))
        return ls


@dataclass
class ResultadoOrbitas:
    energia_final_mev: float
    vueltas: int
    tiempo_total_us: float
    radios_por_vuelta_mm: list
    energias_por_vuelta_mev: list
    f_rev_final_mhz: float
    f_objetivo_mhz: float
    separacion_vueltas_final_mm: float
    deriva_fase_rf_deg: float
    centros_vuelta_mm: list
    dE_turn_medio_kev: float
    motivo_fin: str = "energia_alcanzada"
    trayectoria_xy_mm: list = field(repr=False, default=None)
    estado: str = "simulacion_teorica"
    fuentes: dict = field(default_factory=lambda: FUENTES_ORBITA)
    pendientes: list = field(default_factory=lambda: [
        "Sin mapa OPERA/medido completo: perfil radial = 2 anclas S09 + interpolación CALC",
        "α_dee = 38° es EST declarado; afecta dE/vuelta (±%)",
        "Centrado estático: sin campo de 1er armónico ni bump de centrado",
        "Validación contra órbita medida: BLOQUEADA hasta haz real (memoria 3)",
    ])


def _paso_exacto(campo: CampoKIUBE, x, y, ang, brho, ds):
    B = campo.B(x, y)
    k = B / brho if brho > 0 else 0.0   # convención: giro antihorario (H-)
    dang = k * ds
    if abs(dang) > 1e-13:
        R = 1 / k
        x += R * (math.sin(ang + dang) - math.sin(ang))
        y -= R * (math.cos(ang + dang) - math.cos(ang))
    else:
        x += ds * math.cos(ang)
        y += ds * math.sin(ang)
    return x, y, ang + dang


def simular_aceleracion(campo: CampoKIUBE | None = None,
                        rf: RFKIUBE | None = None,
                        e_inicial_kev: float = 20.0,   # EST declarado
                        e_final_mev: float = 18.0,     # TUP03
                        ds_mm: float = 2.0,
                        fase_inicial_deg: float = 0.0,
                        inicio: str = "equilibrio",
                        max_vueltas: int = 400,
                        max_pasos: int = 6_000_000) -> ResultadoOrbitas:
    """Acelera H- desde la región central hasta e_final con RF real.

    Kicks por potencial de dee real: entrar al dee ΔE = -V_dee(t),
    salir ΔE = +V_dee(t) (q(H-)=-1; ambos dees en fase, h=2 → la separación
    azimutal entrada-salida produce el factor sin(h·α/2) de forma natural).
    inicio="equilibrio": arranca sobre la órbita de equilibrio de E_inicial
    (aproximación estándar; el arranque caótico desde la fuente queda como
    pendiente declarado de la región central).
    Si la fase RF es desfavorable la partícula PIERDE energía: al caer bajo
    e_inicial/4 se declara pérdida por fase (física real de captura central).
    """
    campo = campo or CampoKIUBE()
    rf = rf or RFKIUBE()
    e_piso_kev = e_inicial_kev / 4.0
    motivo_fin = "energia_alcanzada"
    E = e_inicial_kev / 1e3
    if inicio == "equilibrio":
        br0 = brho_mev_tm(E)
        r0 = br0 / campo.b_media(br0 / campo.b0_t * 1e3)  # it. punto fijo
        for _ in range(20):
            r0 = br0 / campo.b_media(r0 * 1e3)
        az0 = math.radians(campo.colina_centro0_deg)
        x, y = r0 * math.cos(az0), r0 * math.sin(az0)
        ang = az0 + math.pi / 2
    else:
        x, y = 2e-3, 0.0                     # fuente (EST)
        ang = math.pi / 2
    t = 0.0
    fase0 = math.radians(fase_inicial_deg)
    gaps = sorted(rf.lineas_gap())
    az_un = math.atan2(y, x)             # azimut desarrollado (sin wrap)
    vueltas = 0
    radios, energias, centros = [], [], []
    tray = [(x * 1e3, y * 1e3)]
    fases_cruce = []
    xs_vuelta, ys_vuelta = [], []
    ds = ds_mm / 1e3

    n_pasos = 0
    while E < e_final_mev and vueltas < max_vueltas and n_pasos < max_pasos:
        n_pasos += 1
        if E * 1e3 < e_piso_kev:
            motivo_fin = "perdida_fase_rf_central"
            break
        ds_eff = min(ds, max(math.hypot(x, y) / 8.0, 2e-4))  # ≤45°/paso en el centro
        v = CLIGHT * math.sqrt(1 - 1 / (1 + E / M_HMINUS) ** 2)
        brho = brho_mev_tm(E)
        x, y, ang = _paso_exacto(campo, x, y, ang, brho, ds_eff)
        t += ds_eff / v
        tray.append((x * 1e3, y * 1e3))
        xs_vuelta.append(x)
        ys_vuelta.append(y)
        az_nuevo = az_un + wrap_pi(math.atan2(y, x) - wrap_pi(az_un))
        # cruce de borde de dee (azimut desarrollado). q(H-) = -1:
        # entrar al dee: ΔE = q·(V_dee - 0); salir: ΔE = q·(0 - V_dee) [keV]
        for g_az, _dee, direccion in gaps:
            g_un = az_un + wrap_pi(g_az - wrap_pi(az_un))
            if az_un < g_un <= az_nuevo:
                v_kv = rf.v_dee_kv * math.cos(rf.omega * t + fase0)
                E += (-direccion) * v_kv / 1e3   # q(H-)=-1: entra → -V, sale → +V
                fases_cruce.append(math.degrees(wrap_pi(rf.omega * t + fase0)))
        az_un = az_nuevo
        if az_un >= math.atan2(tray[0][1], tray[0][0]) + 2 * math.pi * (vueltas + 1):
            vueltas += 1
            radios.append(math.hypot(x, y) * 1e3)
            energias.append(E)
            # centro de la vuelta: ajuste de círculo (Kasa) rápido
            X = np.array(xs_vuelta)
            Y = np.array(ys_vuelta)
            A = np.c_[2 * X, 2 * Y, np.ones(len(X))]
            bvec = X ** 2 + Y ** 2
            try:
                sol, *_ = np.linalg.lstsq(A, bvec, rcond=None)
                centros.append((float(sol[0] * 1e3), float(sol[1] * 1e3)))
            except np.linalg.LinAlgError:
                centros.append((float("nan"), float("nan")))
            xs_vuelta, ys_vuelta = [], []

    if vueltas >= max_vueltas:
        motivo_fin = "max_vueltas"
    elif n_pasos >= max_pasos:
        motivo_fin = "max_pasos"
    if not radios:
        radios, energias = [math.hypot(x, y) * 1e3], [E]
        centros = centros or [(float("nan"), float("nan"))]
    # f_rev final: tiempo de la última vuelta ≈ t_total/vueltas (isócrono)
    t_vuelta = t / max(vueltas, 1)
    f_rev = 1 / t_vuelta
    f_obj = rf.f_mhz * 1e6 / rf.armonico
    dR = radios[-1] - radios[-2] if len(radios) > 1 else float("nan")
    dE_med = (energias[-1] * 1e3 - e_inicial_kev) / max(vueltas, 1)
    fases = np.array(fases_cruce) if fases_cruce else np.zeros(1)
    return ResultadoOrbitas(
        energia_final_mev=energias[-1] if energias else E,
        vueltas=vueltas,
        tiempo_total_us=t * 1e6,
        radios_por_vuelta_mm=[round(r, 3) for r in radios],
        energias_por_vuelta_mev=[round(e, 5) for e in energias],
        f_rev_final_mhz=f_rev / 1e6,
        f_objetivo_mhz=f_obj / 1e6,
        separacion_vueltas_final_mm=dR,
        deriva_fase_rf_deg=float(np.std(fases)) if len(fases) > 1 else 0.0,
        centros_vuelta_mm=[(round(a, 2), round(b, 2)) for a, b in centros],
        dE_turn_medio_kev=dE_med,
        motivo_fin=motivo_fin,
        trayectoria_xy_mm=[(round(a, 2), round(b, 2)) for a, b in tray],
    )


def cargar_perfil_key22(campo: CampoKIUBE, path_csv) -> CampoKIUBE:
    """Carga ⟨B⟩(r) y F(r) MEDIDOS del KIUBE (KEY22, digitalizado) en el campo.

    El CSV (data/kiube/kiube_Bavg_flutter_KEY22.csv) trae huecos: ⟨B⟩ falta en
    el centro (r<23 mm) y en el borde; F falta en 0,285–0,345 m (curvas que se
    cruzan, lectura ambigua — declarado en la cabecera del CSV). Se interpolan
    linealmente sobre los huecos y se declara; en el centro ⟨B⟩ se extiende con
    el primer valor medido (joroba central ~1,033, KEY22 B_centro_rel).
    Además fija gap_valle_mm = 110 mm (KEY22-GEOM: el FE 3D contra el mapa
    medido fijó el valle a 110 mm del plano medio; antes 132 mm estimado).
    """
    d = np.genfromtxt(path_csv, delimiter=",", skip_header=4)
    r_m, b_g, f = d[:, 0], d[:, 1], d[:, 2]
    mb = ~np.isnan(b_g)
    mf = ~np.isnan(f)
    b_full = np.interp(r_m, r_m[mb], b_g[mb])          # huecos interpolados
    f_full = np.interp(r_m, r_m[mf], f[mf])
    campo.key22_r_mm = tuple(float(x * 1e3) for x in r_m)
    campo.key22_bavg_t = tuple(float(x * 1e-4) for x in b_full)  # G → T
    campo.key22_f = tuple(float(x) for x in f_full)
    campo.gap_valle_mm = 110.0
    return campo


def medir_f_rev_orbita(campo: CampoKIUBE, T_mev: float) -> float:
    """f_rev de la órbita cerrada a T_mev integrando un sector (CALC)."""
    brho = brho_mev_tm(T_mev)
    v = CLIGHT * math.sqrt(1 - 1 / (1 + T_mev / M_HMINUS) ** 2)
    c0 = math.radians(campo.colina_centro0_deg)
    r = brho / campo.b0_t
    for _ in range(30):
        r = brho / campo.b_media(r * 1e3)
    x, y = r * math.cos(c0), r * math.sin(c0)
    ang = c0 + math.pi / 2
    t = 0.0
    az = az0 = math.atan2(y, x)
    ds = 5e-5
    for _ in range(4_000_000):
        x, y, ang = _paso_exacto(campo, x, y, ang, brho, ds)
        t += ds / v
        az += wrap_pi(math.atan2(y, x) - wrap_pi(az))
        if az >= az0 + 2 * math.pi / campo.n_sectores:
            break
    return 1.0 / (campo.n_sectores * t)


def calibrar_shim_isocrono(campo: CampoKIUBE, rf: RFKIUBE,
                           T_ref_mev: float = 8.0,
                           iteraciones: int = 2) -> float:
    """Shim multiplicativo de B̄ para f_rev(T_ref) = f_RF/h (CALC declarado).

    Absorbe el efecto de flutter hard-edge sobre el período (la construcción
    B̄=γB0 con campo medio es exacta solo a orden 0 en flutter).
    """
    f_obj = rf.f_mhz * 1e6 / rf.armonico
    for _ in range(iteraciones):
        campo.shim_b0 *= f_obj / medir_f_rev_orbita(campo, T_ref_mev)
    return campo.shim_b0


def energia_equilibrio_mev(campo: CampoKIUBE, r_mm: float) -> float:
    """E de la órbita de equilibrio en r (pc = 0,2998·B̄·r, CALC)."""
    pc = 0.299792458 * campo.b_media(r_mm) * (r_mm / 1e3) * 1e3  # MeV
    return math.sqrt(pc * pc + M_HMINUS ** 2) - M_HMINUS


def calibrar_shim_radial(campo: CampoKIUBE, rf: RFKIUBE,
                         energias_mev=(0.02, 0.1, 0.5, 2.0, 8.0, 18.0),
                         iteraciones: int = 3) -> dict:
    """Tabla de shim shim_r(r) para f_rev(E) = f_RF/h en todo el rango (CALC).

    Cada iteración mide f_rev(E) por órbita cerrada sectorial y multiplica la
    tabla en r_eq(E) por f_obj/f_medida. Converge a |err| < 1e-4 típico.
    """
    f_obj = rf.f_mhz * 1e6 / rf.armonico
    rs = [None] * len(energias_mev)
    for it in range(iteraciones):
        errs = []
        for i, E in enumerate(energias_mev):
            f = medir_f_rev_orbita(campo, E)
            brho = brho_mev_tm(E)
            r = brho / campo.b0_t
            for _ in range(30):
                r = brho / campo.b_media(r * 1e3)
            rs[i] = r * 1e3
            errs.append(f_obj / f)
        if campo.shim_r_mm:
            base = [float(np.interp(r, campo.shim_r_mm, campo.shim_f))
                    for r in rs]
        else:
            base = [1.0] * len(rs)
        campo.shim_r_mm = tuple(rs)
        campo.shim_f = tuple(b * e for b, e in zip(base, errs))
    return {"r_mm": list(campo.shim_r_mm), "shim": list(campo.shim_f),
            "errores_finales": [
                medir_f_rev_orbita(campo, E) / f_obj - 1
                for E in energias_mev]}


def resolver_colina_isocrona(campo: CampoKIUBE | None = None,
                             rf: RFKIUBE | None = None,
                             T_mev: float = 18.0) -> float:
    """Ancho efectivo de colina para f_rev(T) = f_RF/h (órbita cerrada, CALC).

    Mismo criterio que kiube v0.4 (beam.solve_isochronous_hill): integra la
    órbita cerrada en un sector y ajusta el ancho por bisección.
    """
    campo = campo or CampoKIUBE()
    rf = rf or RFKIUBE()
    brho = brho_mev_tm(T_mev)
    f_obj = rf.f_mhz * 1e6 / rf.armonico
    per = 2 * math.pi / campo.n_sectores
    c0 = math.radians(campo.colina_centro0_deg)
    v = CLIGHT * math.sqrt(1 - 1 / (1 + T_mev / M_HMINUS) ** 2)

    def f_rev(ancho_deg: float) -> float:
        campo.colina_ancho_ef_deg = ancho_deg

        def err(r0):
            x, y = r0 * math.cos(c0), r0 * math.sin(c0)
            a = c0 + math.pi / 2
            tt = 0.0
            ds = 2e-4
            for _ in range(200000):
                x, y, a = _paso_exacto(campo, x, y, a, brho, ds)
                tt += ds / v
                if wrap_pi(math.atan2(y, x) - (c0 + per)) >= -1e-6 \
                        and math.atan2(y, x) > c0 + 0.1:
                    break
            return wrap_pi(a - (math.atan2(y, x) + math.pi / 2)), tt

        lo, hi = 0.30, 0.65
        elo, _ = err(lo)
        tq = 0.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            em, tq = err(mid)
            if math.copysign(1, em) == math.copysign(1, elo):
                lo, elo = mid, em
            else:
                hi = mid
        return campo.n_sectores / tq

    lo, hi = 35.0, 65.0
    for _ in range(22):
        m = 0.5 * (lo + hi)
        if f_rev(m) < f_obj:
            lo = m
        else:
            hi = m
    campo.colina_ancho_ef_deg = 0.5 * (lo + hi)
    return campo.colina_ancho_ef_deg
