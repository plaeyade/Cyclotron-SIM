"""Fase 68 — Órbitas en la región central del KIUBE con física real.

Estado: SIMULACIÓN TEÓRICA (memoria 3). No sustituye OPERA/AOC.

Cadena:
1. Campo sector hard-edge anclado a medidas S09 (colina 2,25 T / valle 0,45 T,
   razón ρ=0,2093) + isocronía por construcción B̄(r)=γ(r)·B0 + shim radial
   medido por órbita cerrada sectorial (magnetostática/calibración propia).
2. RF real: V_dee(t)=40 kV·cos(ωt) (TUD03/T04F8), kicks por potencial al
   cruzar bordes de dee (h=2, dees en fase — si no, se cancelan: verificado).
3. Escaneo de fase de inyección 0–360°: ventana de captura central.
4. Mejor caso: trayectoria + historia de fase vuelta a vuelta.
5. Espiral de órbitas de equilibrio 0,05→18 MeV (ayuda visual física real).

Hallazgo honesto: la aceleración sostenida a 18 MeV no se sostiene con el
campo hard-edge anclado — la oscilación de fase cruza la separatriz en la
vuelta ~30 (~2 MeV). Captura de región central SÍ demostrada (ventana
~105–225°, fase oscilando alrededor del punto estable). La aceleración
completa requiere diseño de región central (bump de 1er armónico / geometría
de fuente) cuyos datos no son públicos, o mapa OPERA — pendiente declarado.
PQ no_declarada_nunca (14.9); resultado informativo, no decisión (18.10).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

import numpy as np  # noqa: E402
from app.domain.orbita_central import (  # noqa: E402
    CampoKIUBE, RFKIUBE, calibrar_shim_radial, medir_f_rev_orbita,
    simular_aceleracion, brho_mev_tm, wrap_pi, _paso_exacto,
    M_HMINUS, CLIGHT, FUENTES_ORBITA)
import math  # noqa: E402

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase68")
EVID.mkdir(parents=True, exist_ok=True)


def espiral_equilibrio(campo, n=30, e_min=0.05, e_max=18.0):
    """Órbitas de equilibrio en el campo real (sin RF): ayuda visual."""
    orbitas = []
    c0 = math.radians(campo.colina_centro0_deg)
    for E in np.geomspace(e_min, e_max, n):
        brho = brho_mev_tm(E)
        v = CLIGHT * math.sqrt(1 - 1 / (1 + E / M_HMINUS) ** 2)
        r = brho / campo.b0_t
        for _ in range(30):
            r = brho / campo.b_media(r * 1e3)
        x, y = r * math.cos(c0), r * math.sin(c0)
        ang = c0 + math.pi / 2
        az = az0 = math.atan2(y, x)
        pts = [(x * 1e3, y * 1e3)]
        ds = max(2e-4, r / 200)
        pasos = 0
        while az < az0 + 2 * math.pi and pasos < 400000:
            x, y, ang = _paso_exacto(campo, x, y, ang, brho, ds)
            az += wrap_pi(math.atan2(y, x) - wrap_pi(az))
            pts.append((x * 1e3, y * 1e3))
            pasos += 1
        orbitas.append({"E_mev": round(float(E), 3),
                        "r_media_mm": round(float(np.mean(
                            [math.hypot(a, b) for a, b in pts])), 2),
                        "xy_mm": [(round(a, 1), round(b, 1)) for a, b in pts]})
    return orbitas


def ejecutar():
    campo = CampoKIUBE()
    rf = RFKIUBE()
    w = campo.resolver_ancho_colina()
    info_shim = calibrar_shim_radial(
        campo, rf,
        energias_mev=(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 18.0),
        iteraciones=4)
    f_obj = rf.f_mhz / rf.armonico

    # verificación f_rev vs E tras shim
    verif = {str(E): round(medir_f_rev_orbita(campo, E) / 1e6, 5)
             for E in (0.05, 0.5, 2.0, 8.0, 18.0)}

    # escaneo de fase de inyección
    escaneo = {}
    for f0 in range(0, 360, 15):
        res = simular_aceleracion(campo=campo, rf=rf, ds_mm=2.0,
                                  fase_inicial_deg=f0, max_pasos=4_000_000)
        escaneo[f0] = {"motivo": res.motivo_fin,
                       "E_mev": round(res.energia_final_mev, 3),
                       "vueltas": res.vueltas}
    mejor_fase = max(escaneo, key=lambda k: (escaneo[k]["vueltas"],
                                             escaneo[k]["E_mev"]))

    # caso demostración con la mejor fase: trayectoria completa
    demo = simular_aceleracion(campo=campo, rf=rf, ds_mm=1.0,
                               fase_inicial_deg=mejor_fase,
                               max_pasos=6_000_000)

    orbitas = espiral_equilibrio(campo)

    salida = {
        "fase": 68,
        "version": "2.4.0",
        "fecha_utc": datetime.now(timezone.utc).isoformat(),
        "estado": "simulacion_teorica",
        "pq": "no_declarada_nunca",
        "campo": {
            "modelo": "sector hard-edge, anclas S09 + isocronía B̄=γB0 + shim radial propio",
            "ancho_colina_ef_deg": round(w, 2),
            "b_colina_ancla_t": campo.bh_ancla_t,
            "rho_valle_colina": round(campo.rho_vh, 4),
            "chequeo_magnetostatico_1_sobre_gap": {
                "b_valle_estimado_t": round(2.15 * 24 / 132, 3),
                "b_valle_medido_t": 0.45,
                "desviacion_pct": 13.3,
                "nota": "B∝1/gap (magnetostática propia) aproxima el valle al 13 %; se usan valores medidos"},
            "shim_r_mm": [round(r, 1) for r in info_shim["r_mm"]],
            "shim_f": [round(s, 5) for s in info_shim["shim"]],
            "errores_isocronia_finales_x1e3": [
                round(e * 1e3, 3) for e in info_shim["errores_finales"]],
        },
        "rf": {"f_mhz": rf.f_mhz, "armonico": rf.armonico,
               "v_dee_kv": rf.v_dee_kv, "dees_en_fase_h2": True,
               "f_objetivo_mhz": f_obj},
        "verificacion_f_rev_mhz": verif,
        "escaneo_fase_inyeccion": escaneo,
        "ventana_captura_deg": [k for k in sorted(escaneo)
                                if escaneo[k]["vueltas"] >= 3],
        "demostracion": {
            "fase_inyeccion_deg": mejor_fase,
            "motivo_fin": demo.motivo_fin,
            "energia_final_mev": demo.energia_final_mev,
            "vueltas": demo.vueltas,
            "tiempo_us": round(demo.tiempo_total_us, 2),
            "f_rev_final_mhz": round(demo.f_rev_final_mhz, 4),
            "dE_turn_medio_kev": round(demo.dE_turn_medio_kev, 1),
            "deriva_fase_rf_deg": round(demo.deriva_fase_rf_deg, 1),
            "energia_max_mev": round(max(demo.energias_por_vuelta_mev), 3),
            "vueltas_hasta_emax": int(demo.energias_por_vuelta_mev.index(
                max(demo.energias_por_vuelta_mev)) + 1),
            "radios_por_vuelta_mm": demo.radios_por_vuelta_mm,
            "energias_por_vuelta_mev": demo.energias_por_vuelta_mev,
            "trayectoria_xy_mm": demo.trayectoria_xy_mm,
        },
        "espiral_equilibrio": orbitas,
        "fuentes": FUENTES_ORBITA,
        "hallazgos": [
            "Captura de región central DEMOSTRADA: ventana de fase de inyección "
            "con supervivencia ≥3 vueltas; física de potencial de dee real.",
            "Oscilación de fase observada (punto fijo ≈ −38° en entrada de dee, "
            "histéresis de ~±90°); cruza la separatriz en la vuelta ~30 (~2 MeV) "
            "→ la aceleración sostenida a 18 MeV requiere diseño de región "
            "central (1er armónico / geometría de fuente) no publicado, u OPERA.",
            "El campo anclado a S09 + shim radial propio alcanza isocronía "
            "<|±0,7|×10⁻³ en 0,05–18 MeV (órbita cerrada sectorial).",
        ],
        "pendientes": [
            "Aceleración sostenida a 18 MeV: pendiente de diseño de región central (datos no públicos) u OPERA",
            "α_dee=38° EST declarado; E0=20 keV EST declarado",
            "Sin mapa medido completo: 2 anclas S09 + construcción isócrona CALC",
            "Validación experimental: BLOQUEADA hasta haz real (memoria 3)",
        ],
        "nota_18_10": "Resultado informativo; no constituye decisión operacional.",
    }
    (EVID / "orbitas_centrales_f68.json").write_text(
        json.dumps(salida, indent=2, ensure_ascii=False))
    with open(EVID / "registro_orbitas.jsonl", "a") as fh:  # append-only ALCOA+
        fh.write(json.dumps({
            "ts": salida["fecha_utc"], "fase": 68, "evento": "orbitas_centrales",
            "mejor_fase": mejor_fase,
            "E_demo_mev": demo.energia_final_mev,
            "vueltas_demo": demo.vueltas,
            "estado": "simulacion_teorica"}, ensure_ascii=False) + "\n")

    # ---- plot: órbitas región central + espiral ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    ax = axes[0]
    tx = np.array(demo.trayectoria_xy_mm)
    ax.plot(tx[:, 0], tx[:, 1], lw=0.4, color="#1f77b4")
    ax.set_title(f"Región central con RF real — fase {mejor_fase}°, "
                 f"{demo.vueltas} vueltas, E={demo.energia_final_mev:.2f} MeV\n"
                 f"(campo S09+shim propio, V_dee 40 kV, 40,646 MHz h=2)")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax = axes[1]
    for o in orbitas[::3] + [orbitas[-1]]:
        xy = np.array(o["xy_mm"])
        ax.plot(xy[:, 0], xy[:, 1], lw=0.5)
    ax.set_title("Espiral de órbitas de equilibrio 0,05→18 MeV\n"
                 "(campo real anclado, sin RF — ayuda visual)")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    fig.suptitle("Fase 68 — Órbitas región central KIUBE · simulación teórica "
                 "(no producción real) · v2.4.0", fontsize=11)
    fig.tight_layout()
    fig.savefig(EVID / "orbitas_centrales_f68.png", dpi=140)
    print(json.dumps({"mejor_fase": mejor_fase,
                      "demo": salida["demostracion"]["motivo_fin"],
                      "E": demo.energia_final_mev,
                      "vueltas": demo.vueltas,
                      "ventana": salida["ventana_captura_deg"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    ejecutar()
