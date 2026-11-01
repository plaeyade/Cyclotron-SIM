"""
Pérdida de haz H⁻ por stripping en gas residual (Escalón A — F65).

Cadena física: P(r) declarada → n_gas(r) → σ_stripping(E) citada →
pérdida por vuelta a lo largo de la órbita → H⁰ sale en trayectoria
tangente rectilínea → impacto en la pared de la cámara (mapa azimutal).

Fuentes de las secciones eficaces (regla 49.11, cita verbatim):
  - PLUM-CYR: M. Plum, "Beam Loss in Linacs", CERN Yellow Reports
    (e-publishing.cern.ch, art. 229), ec. (1)-(2):
        σ_N2/O2 = 7×10⁻¹⁹/β² cm² por átomo
        σ_H2    = 1×10⁻¹⁹/β² cm² por átomo
    basadas en G. Gillespie, Phys. Rev. A 15 (1977) 563;
    Phys. Rev. A 16 (1977) 943. Stripping doble ≈ 4 % del simple.
  - NAKHAI2022: Nakhai-Pour et al., "Cyclotron Vacuum Model and H⁻ Gas
    Stripping Loss" (2022): composición de gas residual en ciclotrón
    70 % H₂O, 30 % aire, 0.15 % gas raro.
  - KIUBE-V04: params.yaml del modelo kiube v0.4 (ΔE = 98,5 keV/vuelta,
    ~183 vueltas, <B>(18 MeV) = 1,360 T, cámara r = 548 mm) — medidas
    T04F2/CALC declaradas por el propio modelo.

σ_H2O no está tabulada en PLUM-CYR: se deriva por aditividad atómica
σ_H2O ≈ σ_O + 2·σ_H (CALC declarado: O = 7×10⁻¹⁹/β², H = 0,5×10⁻¹⁹/β²
por átomo). Es la aproximación estándar de mezcla atómica; su
incertidumbre queda subsumida en la incertidumbre de P(r).

Estado canónico de todo resultado: `simulacion_teorica`. P(r) es un
perfil parametrizado DECLARADO (supuesto_declarado_pendiente_fuente);
será sustituido por el mapa Molflow+ en el Escalón B. Nada aquí habilita
producción clínica (18.19 pendiente_entorno; PQ no_declarada_nunca 14.9).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# --- constantes físicas ---
K_B = 1.380649e-23          # J/K
M_H_MEV = 938.272           # masa protón (H⁻ ≈ p), MeV/c²
MBAR_A_PA = 100.0           # 1 mbar = 100 Pa (F82: corregido; hasta v3.7.0 valía 0,1 → densidad y pérdidas 1000× bajas)

# --- secciones eficaces citadas: sigma = C / beta² (cm² por átomo) ---
_C_ATOMO = {  # coeficiente C en 1e-19 cm²
    "H": (0.5e-19, "PLUM-CYR ec.(2): H2 = 1e-19/beta² cm²/átomo-molécula → 0,5e-19 por átomo H"),
    "N": (7.0e-19, "PLUM-CYR ec.(1): N2/O2 = 7e-19/beta² cm²/átomo"),
    "O": (7.0e-19, "PLUM-CYR ec.(1): N2/O2 = 7e-19/beta² cm²/átomo"),
}

# moléculas como suma de átomos (aditividad, CALC declarado)
_MOLECULA_ATOMOS = {
    "H2":  {"H": 2},
    "N2":  {"N": 2},
    "O2":  {"O": 2},
    "H2O": {"O": 1, "H": 2},   # CALC: aditividad atómica sobre PLUM-CYR
    "aire": {"N": 1.56, "O": 0.42},  # 78% N2 + 21% O2 por molécula "media"
}

COMPOSICION_RESIDUAL_NAKHAI = {"H2O": 0.70, "aire": 0.30}  # NAKHAI2022 (0.15% gases raros omitidos, renormalizado)

FUENTES = {
    "secciones_stripping": "PLUM-CYR: M. Plum, 'Beam Loss in Linacs', CERN Yellow Reports (art. 229), ec. (1)-(2); base: Gillespie, Phys. Rev. A 15 (1977) 563 y 16 (1977) 943",
    "composicion_gas": "NAKHAI2022: Nakhai-Pour et al., 'Cyclotron Vacuum Model and H- Gas Stripping Loss' (2022): 70% H2O, 30% aire",
    "orbita_kiube": "KIUBE-V04: params.yaml kiube v0.4 (ΔE=98,5 keV/vuelta CALC, <B>=1,360 T CALC, cámara r=548 mm T04F2)",
    "h2o_aditividad": "CALC declarado: σ_H2O = σ_O + 2·σ_H por aditividad atómica",
}

FACTOR_STRIPPING_DOBLE = 0.04  # PLUM-CYR: doble ≈ 4 % del simple (H⁻→H⁺)


def beta2(energia_mev: float) -> float:
    """β² relativista para H⁻ (≈ protón)."""
    g = 1.0 + energia_mev / M_H_MEV
    return 1.0 - 1.0 / (g * g)


def sigma_molecula_cm2(gas: str, energia_mev: float) -> float:
    """σ de stripping simple H⁻→H⁰ para una molécula de gas (cm²)."""
    b2 = beta2(energia_mev)
    atomos = _MOLECULA_ATOMOS[gas]
    return sum(_C_ATOMO[a][0] * n for a, n in atomos.items()) / b2


def sigma_mezcla_cm2(composicion: dict[str, float], energia_mev: float) -> float:
    """σ efectiva de la mezcla residual (cm² por molécula)."""
    return sum(f * sigma_molecula_cm2(g, energia_mev) for g, f in composicion.items())


def densidad_gas_cm3(presion_mbar: float, temperatura_k: float = 303.0) -> float:
    """n = P/(kB·T) en moléculas/cm³. T = 303 K (ejemplo PLUM-CYR)."""
    return presion_mbar * MBAR_A_PA / (K_B * temperatura_k) / 1e6


def radio_orbita_mm(energia_mev: float, b_media_t: float) -> float:
    """r = p/(qB) en mm, p de T relativista. KIUBE-V04: <B>(18 MeV)=1,360 T."""
    p_mev_c = math.sqrt(energia_mev * (energia_mev + 2 * M_H_MEV))
    # r[m] = p[GeV/c]/(0.2998·B[T]); con p en MeV/c queda directamente en mm
    return p_mev_c / (0.299792458 * b_media_t)


@dataclass
class PerfilPresion:
    """P(r) parametrizado DECLARADO o desde mapa TPMC (vacuum_tpmc, F66).

    modelos:
      - "uniforme": P(r) = p_base
      - "radial":   P(r) = p_base + (p_centro - p_base)·(1 - r/r_pared)
        (la fuente de H₂ de las fuentes de iones está en el centro;
        las 4 bombas ODP evacúan en la periferia → máximo en el centro)
      - "mapa":     interpolación lineal sobre (mapa_r_mm, mapa_p_mbar)
        proveniente del motor TPMC (F66) — mapa 3D real, no parametrizado
    """
    p_base_mbar: float
    modelo: str = "uniforme"
    p_centro_mbar: float | None = None
    r_pared_mm: float = 548.0
    mapa_r_mm: tuple | None = None
    mapa_p_mbar: tuple | None = None

    @classmethod
    def desde_mapa(cls, r_mm, p_mbar):
        return cls(p_base_mbar=p_mbar[-1], modelo="mapa",
                   mapa_r_mm=tuple(r_mm), mapa_p_mbar=tuple(p_mbar))

    def presion_mbar(self, r_mm: float) -> float:
        if self.modelo == "mapa" and self.mapa_r_mm:
            import bisect
            rs, ps = self.mapa_r_mm, self.mapa_p_mbar
            if r_mm <= rs[0]:
                return ps[0]
            if r_mm >= rs[-1]:
                return ps[-1]
            i = bisect.bisect_right(rs, r_mm) - 1
            t = (r_mm - rs[i]) / (rs[i + 1] - rs[i])
            return ps[i] * (1 - t) + ps[i + 1] * t
        if self.modelo == "uniforme" or self.p_centro_mbar is None:
            return self.p_base_mbar
        f = max(0.0, 1.0 - r_mm / self.r_pared_mm)
        return self.p_base_mbar + (self.p_centro_mbar - self.p_base_mbar) * f


@dataclass
class ResultadoStripping:
    fraccion_perdida: float
    perdida_por_vuelta: list[float]
    impactos_azimut: list[float]     # histograma 36 bins (10°), fracción
    potencia_perdida_w: float
    n_vueltas: int
    parametros: dict = field(default_factory=dict)
    fuentes: dict = field(default_factory=lambda: dict(FUENTES))
    estado: str = "simulacion_teorica"
    pendientes: tuple = (
        "P(r) parametrizado: sustituir por mapa Molflow+ (Escalón B)",
        "σ(E) por gas y energía exacta: tablas Gillespie completas",
        "geometría 2D plano medio: impactos z reales con modelo 3D",
    )


def simular_perdida_gas(
    corriente_ua: float,
    energia_final_mev: float,
    delta_e_kev_vuelta: float,
    b_media_t: float,
    perfil: PerfilPresion,
    composicion: dict[str, float] = COMPOSICION_RESIDUAL_NAKHAI,
    r_pared_mm: float = 548.0,
    incluir_doble: bool = True,
) -> ResultadoStripping:
    """Integra la pérdida vuelta a vuelta y traza los H⁰ a la pared.

    Órbita: espiral de r(E) = p/(qB) (KIUBE-V04). Pérdida por vuelta i:
    ΔI/I = σ_mezcla(E_i)·n(P(r_i))·2πr_i. El H⁰ sale tangente desde
    (r_i, φ_i) con φ_i uniforme en [0, 2π) por vuelta (espiral densa) y
    choca con el cilindro de la cámara en r = r_pared.
    """
    n_vueltas = max(1, round(energia_final_mev * 1e3 / delta_e_kev_vuelta))
    supervivencia = 1.0
    perdida_vuelta: list[float] = []
    impactos = [0.0] * 36
    energia_perdida_mev = 0.0
    factor = 1.0 + (FACTOR_STRIPPING_DOBLE if incluir_doble else 0.0)

    for i in range(1, n_vueltas + 1):
        e = energia_final_mev * i / n_vueltas
        r_mm = radio_orbita_mm(e, b_media_t)
        if r_mm >= r_pared_mm:
            r_mm = r_pared_mm * 0.999  # última vuelta junto a la pared
        p_mbar = perfil.presion_mbar(r_mm)
        n_cm3 = densidad_gas_cm3(p_mbar)
        sig = sigma_mezcla_cm2(composicion, e)
        camino_cm = 2.0 * math.pi * (r_mm / 10.0)
        frac = min(sig * n_cm3 * camino_cm * factor, 1.0)
        perdida = supervivencia * frac
        supervivencia -= perdida
        perdida_vuelta.append(perdida)
        energia_perdida_mev += perdida * e
        # trazado H⁰: punto de stripping uniforme en la vuelta; el neutro
        # sigue la tangente hasta el cilindro. Impacto en azimut φ_imp.
        phi_s = 2.0 * math.pi * (i * 0.6180339887 % 1.0)  # muestreo áureo
        # geometría: cuerda tangente desde (r, φ_s) al círculo R:
        dphi = math.acos(min(r_mm / r_pared_mm, 1.0))
        phi_imp = (phi_s + dphi) % (2.0 * math.pi)
        impactos[int(phi_imp / (2.0 * math.pi) * 36) % 36] += perdida

    fraccion_total = 1.0 - supervivencia
    # potencia: corriente incidente × energía media ponderada de los H⁰
    # perdidos (fracción·E ya acumulada): P = I·<E>_ponderado
    potencia_w = corriente_ua * 1e-6 * (energia_perdida_mev * 1e6)  # I·V = W

    total_imp = sum(impactos) or 1.0
    return ResultadoStripping(
        fraccion_perdida=fraccion_total,
        perdida_por_vuelta=perdida_vuelta,
        impactos_azimut=[x / total_imp for x in impactos],
        potencia_perdida_w=potencia_w,
        n_vueltas=n_vueltas,
        parametros={
            "corriente_ua": corriente_ua,
            "energia_final_mev": energia_final_mev,
            "delta_e_kev_vuelta": delta_e_kev_vuelta,
            "b_media_t": b_media_t,
            "r_pared_mm": r_pared_mm,
            "presion": {"modelo": perfil.modelo, "p_base_mbar": perfil.p_base_mbar,
                        "p_centro_mbar": perfil.p_centro_mbar},
            "composicion": composicion,
            "stripping_doble_incluido": incluir_doble,
            "energia_media_h0_mev": (energia_perdida_mev / fraccion_total) if fraccion_total else 0.0,
        },
    )
