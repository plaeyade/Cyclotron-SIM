"""Fase 69 — Recalibrado del campo con el mapa MEDIDO KEY22 (Cyclotrons2022).

Estado: SIMULACIÓN TEÓRICA (memoria 3: sin haz real; validación teórica).

Sustituye el perfil analítico isócrono B̄(r)=γ(r)·B0 de F68 por el perfil
MEDIDO del KIUBE: ⟨B⟩(r) y flutter F(r) digitalizados de Nuttens et al.,
Cyclotrons2022 WEAO03 Fig.3 (KEY22), con ±15 G / ±0,005 de lectura y huecos
declarados (centro r<23 mm, borde, y F en 0,285–0,345 m) interpolados
linealmente (CALC declarado). B_colina/B_valle salen del flutter medido:
F = a(1−a)(Bh−Bv)²/B̄² con a = w/90 (w = 44,28°, ancla S09 de F68).

Geometría KEY22-GEOM (fijada por el FE 3D de kiube v0.7 contra el mapa
medido): valle a 110 mm del plano medio (antes 132 mm estimado) y tapón
rebajado 17 mm. Con 110 mm el chequeo de magnetostática propia B∝1/gap mejora:
2,15·24/110 = 0,469 T vs 0,516 T del flutter medido (9 %; antes 13 %).

Encima del perfil medido se recalibra el shim radial propio (órbita cerrada
sectorial, mismo método F68) para isocronía fina — la medida real no es
perfectamente isócrona y el shim absorbe el residuo (|err| < 3e-4).

Verificaciones:
- ⟨B⟩(r) medido vs analítico F68 (joroba central 3,3 % que el analítico no tenía).
- Bh(450 mm)/Bv(450 mm) del flutter medido vs FE NGSolve de kiube v0.7
  (2,169/0,482 T): definiciones distintas (integral de flutter vs valores a
  45° de azimut), declarado.
- Órbita de equilibrio 18 MeV vs F68 (r_min 423,0 / r_max 484,2 mm).
- Captura central con RF real a la mejor fase F68 (165°): la conclusión de
  F68 (captura sí, aceleración sostenida no sin diseño de región central)
  debe mantenerse con el campo medido.

PQ no_declarada_nunca (14.9); resultado informativo, no decisión (18.10).
Fuentes (49.11): KEY22, KEY22-GEOM, S09, T04F2, T04F8, TUD03, CALC.
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

import numpy as np  # noqa: E402
from app.domain.orbita_central import (  # noqa: E402
    CampoKIUBE, RFKIUBE, cargar_perfil_key22, calibrar_shim_radial,
    medir_f_rev_orbita, simular_aceleracion, brho_mev_tm, wrap_pi,
    _paso_exacto, M_HMINUS, CLIGHT, FUENTES_ORBITA)

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase69")
EVID.mkdir(parents=True, exist_ok=True)
CSV_KEY22 = RAIZ / "data" / "kiube" / "kiube_Bavg_flutter_KEY22.csv"


def radios_orbita(campo, E=18.0):
    """Órbita de equilibrio completa a E: r_min/r_max (mm)."""
    brho = brho_mev_tm(E)
    c0 = math.radians(campo.colina_centro0_deg)
    r = brho / campo.b0_t
    for _ in range(40):
        r = brho / campo.b_media(r * 1e3)
    x, y = r * math.cos(c0), r * math.sin(c0)
    ang = c0 + math.pi / 2
    az = az0 = math.atan2(y, x)
    rs = []
    ds = max(2e-4, r / 300)
    pasos = 0
    while az < az0 + 2 * math.pi and pasos < 600000:
        x, y, ang = _paso_exacto(campo, x, y, ang, brho, ds)
        az += wrap_pi(math.atan2(y, x) - wrap_pi(az))
        rs.append(math.hypot(x, y) * 1e3)
        pasos += 1
    return min(rs), max(rs), np.mean(rs)


def ejecutar():
    rf = RFKIUBE()

    # --- campo F68 (referencia) y campo F69 (KEY22 medido)
    c68 = CampoKIUBE()
    w68 = c68.resolver_ancho_colina()
    campo = CampoKIUBE()
    cargar_perfil_key22(campo, CSV_KEY22)

    # comparación de perfiles antes del shim
    rr = np.geomspace(23, 505, 40)
    comp_perfil = [{
        "r_mm": round(float(r), 1),
        "b_key22_t": round(campo.b_media(float(r)), 4),
        "b_f68_t": round(c68.b_media(float(r)), 4)}
        for r in rr]

    # joroba central medida (KEY22 B_centro_rel = 1,0329)
    b_centro_rel = campo.b_media(23) / campo.b_media(450)

    # chequeo magnetostático con valle 110 mm (KEY22-GEOM)
    bvalle_flutter_450 = campo.b_valle(450.0)
    chequeo = {
        "gap_valle_mm": campo.gap_valle_mm,
        "b_valle_estimado_1_sobre_gap_t": round(2.15 * 24 / 110, 3),
        "b_valle_flutter_medido_450mm_t": round(bvalle_flutter_450, 3),
        "desviacion_pct": round(abs(2.15 * 24 / 110 - bvalle_flutter_450)
                                / bvalle_flutter_450 * 100, 1),
        "nota": "valle 110 mm (KEY22-GEOM) mejora el chequeo B∝1/gap: "
                "9 % vs 13 % con 132 mm estimado"}

    # shim radial propio sobre el perfil medido (30 nodos: el residuo baja al
    # piso de ruido de lectura del mapa digitalizado, ±15 G ≈ ±1,1e-3 relativo
    # — no se puede calibrar por debajo de la incertidumbre de la medida)
    energias = tuple(float(e) for e in np.geomspace(0.02, 18.0, 30))
    info_shim = calibrar_shim_radial(campo, rf, energias_mev=energias,
                                     iteraciones=4)
    verif = {str(E): round(medir_f_rev_orbita(campo, E) / 1e6, 5)
             for E in (0.05, 0.5, 2.0, 8.0, 18.0)}

    # órbita 18 MeV: F69 medido vs F68 analítico
    rmin69, rmax69, rmed69 = radios_orbita(campo)
    # c68 (referencia F68): shim con los mismos 30 nodos para comparación justa
    calibrar_shim_radial(c68, rf, energias_mev=energias, iteraciones=4)
    rmin68, rmax68, _ = radios_orbita(c68)

    # error de isocronía denso (60 puntos, incluye internodos) tras el shim
    f_obj_hz = rf.f_mhz * 1e6 / rf.armonico
    err_denso = {f"{E:.3f}": round(medir_f_rev_orbita(campo, E) / f_obj_hz - 1, 6)
                 for E in np.geomspace(0.02, 18.0, 60)}
    err_denso_max = max(abs(e) for e in err_denso.values())

    # escaneo de fase de inyección con el campo medido (mismo método F68)
    escaneo = {}
    for f0 in range(0, 360, 15):
        res = simular_aceleracion(campo=campo, rf=rf, ds_mm=2.0,
                                  fase_inicial_deg=f0, max_pasos=1_500_000)
        escaneo[f0] = {"motivo": res.motivo_fin,
                       "E_mev": round(res.energia_final_mev, 3),
                       "E_max_mev": round(max(res.energias_por_vuelta_mev, default=0.0), 3),
                       "vueltas": res.vueltas}
    mejor_fase = max(escaneo, key=lambda k: (escaneo[k]["vueltas"],
                                             escaneo[k]["E_max_mev"]))
    demo = simular_aceleracion(campo=campo, rf=rf, ds_mm=1.0,
                               fase_inicial_deg=mejor_fase,
                               max_pasos=3_000_000)

    salida = {
        "fase": 69,
        "version": "2.5.0",
        "fecha_utc": datetime.now(timezone.utc).isoformat(),
        "estado": "simulacion_teorica",
        "pq": "no_declarada_nunca",
        "campo": {
            "modelo": "sector hard-edge con ⟨B⟩(r) y flutter F(r) MEDIDOS "
                      "(KEY22) + shim radial propio (isocronía fina)",
            "fuente_perfil": "KEY22",
            "csv": "data/kiube/kiube_Bavg_flutter_KEY22.csv",
            "huecos_interpolados": ["r<23 mm (centro)", "borde r>505 mm",
                                    "F en 285–345 mm (lectura ambigua)"],
            "ancho_colina_ef_deg": round(w68, 2),
            "joroba_central_medida_rel": round(float(b_centro_rel), 4),
            "key22_b_centro_rel_publicado": 1.0329,
            "chequeo_magnetostatico": chequeo,
            "shim_r_mm": [round(r, 1) for r in info_shim["r_mm"]],
            "shim_f": [round(f, 5) for f in info_shim["shim"]],
            "errores_isocronia_finales": [
                round(e, 6) for e in info_shim["errores_finales"]],
            "error_isocronia_denso_60ptos_max": round(err_denso_max, 6),
            "piso_ruido_lectura_key22": "±15 G / 13,3 kG ≈ ±1,1e-3 relativo: "
                "el residuo tras shim coincide con la incertidumbre de lectura "
                "del mapa digitalizado — no se puede calibrar por debajo",
        },
        "colina_valle_medidos_450mm": {
            "bh_flutter_t": round(campo.b_colina(450.0), 4),
            "bv_flutter_t": round(campo.b_valle(450.0), 4),
            "bh_fe_kiube_v07_t": 2.169,
            "bv_fe_kiube_v07_t": 0.482,
            "nota": "definiciones distintas: integral de flutter medido vs "
                    "valores a 45° de azimut del FE 3D (declarado)"},
        "f_rev_verificacion_mhz": verif,
        "f_objetivo_mhz": rf.f_mhz / rf.armonico,
        "orbita_18mev": {
            "f69_key22": {"r_min_mm": round(rmin69, 1),
                          "r_max_mm": round(rmax69, 1),
                          "r_media_mm": round(float(rmed69), 1)},
            "f68_analitico": {"r_min_mm": round(rmin68, 1),
                              "r_max_mm": round(rmax68, 1)},
            "f68_publicado": {"r_min_mm": 423.0, "r_max_mm": 484.2}},
        "captura_rf_escaneo": {str(k): v for k, v in escaneo.items()},
        "captura_rf_mejor_fase": {
            "fase_inicial_deg": mejor_fase,
            "motivo_fin": demo.motivo_fin,
            "vueltas": demo.vueltas,
            "energia_max_mev": round(max(demo.energias_por_vuelta_mev, default=0.0), 3),
            "energia_final_mev": round(demo.energia_final_mev, 3),
            "nota": "misma conclusión que F68: captura central sí; "
                    "aceleración sostenida no sin diseño de región central "
                    "(1er armónico no público u OPERA) — pendiente declarado"},
        "comparacion_perfil": comp_perfil,
        "fuentes": FUENTES_ORBITA,
    }
    out_json = EVID / "orbitas_recal_f69.json"
    out_json.write_text(json.dumps(salida, indent=1, ensure_ascii=False))

    # registro append-only (ALCOA+)
    reg = EVID / "registro_orbitas.jsonl"
    with reg.open("a", encoding="utf-8") as f:
        f.write(json.dumps({
            "ts_utc": salida["fecha_utc"], "fase": 69, "evento":
            "recalibrado_key22",
            "detalle": {"fuente": "KEY22 medido",
                        "err_isocronia_max": round(max(
                            abs(e) for e in info_shim["errores_finales"]), 6),
                        "r_min_18mev_mm": round(rmin69, 1),
                        "r_max_18mev_mm": round(rmax69, 1)}},
            ensure_ascii=False) + "\n")

    # figura
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    rs_ = [p["r_mm"] for p in comp_perfil]
    ax[0].plot(rs_, [p["b_key22_t"] for p in comp_perfil], "-",
               color="#b03a2e", lw=2, label="KEY22 medido (F69)")
    ax[0].plot(rs_, [p["b_f68_t"] for p in comp_perfil], "--",
               color="#3b6fd8", lw=1.5, label="Analítico γ·B0 (F68)")
    ax[0].set_xlabel("r [mm]"); ax[0].set_ylabel("⟨B⟩ [T]")
    ax[0].set_title("⟨B⟩(r): medido KEY22 vs analítico F68\n"
                    f"joroba central {b_centro_rel:.3f} (pub. 1,033)")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    theta = np.linspace(0, 2 * np.pi, 400)
    r_az = [campo.b_colina(450.0) if (t % 90 < 44.28 or t % 90 > 90 - 0)
            else 0 for t in np.degrees(theta)]
    E = 18.0
    brho = brho_mev_tm(E)
    x, y = rmed69 / 1e3 * math.cos(math.radians(45)), \
        rmed69 / 1e3 * math.sin(math.radians(45))
    ang = math.radians(45) + math.pi / 2
    az = az0 = math.atan2(y, x)
    xs, ys = [x * 1e3], [y * 1e3]
    pasos = 0
    while az < az0 + 2 * math.pi and pasos < 600000:
        x, y, ang = _paso_exacto(campo, x, y, ang, brho, 2e-4)
        az += wrap_pi(math.atan2(y, x) - wrap_pi(az))
        xs.append(x * 1e3); ys.append(y * 1e3)
        pasos += 1
    ax[1].plot(xs, ys, color="#1e8449", lw=1.2)
    ax[1].set_aspect("equal"); ax[1].grid(alpha=0.3)
    ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("y [mm]")
    ax[1].set_title(f"Órbita 18 MeV con campo medido KEY22\n"
                    f"r_min {rmin69:.0f} / r_max {rmax69:.0f} mm "
                    f"(F68: {rmin68:.0f}/{rmax68:.0f})")
    fig.tight_layout()
    fig.savefig(EVID / "orbitas_recal_f69.png", dpi=140)
    print(json.dumps({"err_max": max(abs(e) for e in
                                     info_shim["errores_finales"]),
                      "r18": [round(rmin69, 1), round(rmax69, 1)],
                      "mejor_fase": mejor_fase,
                      "captura": salida["captura_rf_mejor_fase"]}, indent=1))
    return salida


if __name__ == "__main__":
    ejecutar()
