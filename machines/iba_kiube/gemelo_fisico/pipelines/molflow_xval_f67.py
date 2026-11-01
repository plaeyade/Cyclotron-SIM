"""Fase 67 — Validación cruzada Molflow oficial (CERN) vs TPMC propio (F66).

Estado: SIMULACIÓN TEÓRICA (memoria 3: sin haz real; validación teórica).
Cierre del pendiente declarado en F66: cotejar el mapa de presión propio
contra el motor oficial Molflow+ 2.11.1 (binario molflowCLI, CERN).

Metodología:
1. tools/kiube_camara_molflow.py genera kiube_camara.xml (cámara KIUBE-V04,
   R=548 mm, gap=132 mm, 4 puertos r=127 mm @375 mm, sticking=0.25,
   fuente H2 central + desgasificación de pared 1e-10 mbar·l/s/cm²).
2. molflowCLI -d 2e6 --writeFacetQuantities -> facet_physics.csv.
3. Binning radial de presión de facetas de tapas (|z|>60 mm).
4. Comparación con mapa TPMC propio (fase66/mapa_vacio_f66.json).
5. Pérdida por stripping con mapa Molflow vs mapa propio.

Fuentes (49.11):
- MOLFLOW-2.11.1: binario oficial molflow_fedora_2.11.1.zip (usuario, 2026-10).
- ADY-IPAC14: Ady & Kersevan, IPAC'14 WEPME038 (algoritmo TPMC).
- KIUBE-V04: parámetros geométricos declarados en F65/F66.
Figuras sin fuente -> supuesto_declarado_pendiente_fuente (Informe 49.11).
PQ: no_declarada_nunca (Informe 14.9). Recomendación != decisión (18.10).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))  # app/domain vive en gemelo_fisico/app (antes services/yield_engine de RPOS)
for _m in [k for k in list(sys.modules) if k == "app" or k.startswith("app.")]:
    del sys.modules[_m]

from app.domain.gas_stripping import PerfilPresion, simular_perdida_gas  # noqa: E402

EVID = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase67")
EVID.mkdir(parents=True, exist_ok=True)

FUENTES = {
    "MOLFLOW-2.11.1": "Binario oficial molflowCLI 2.11.1 (CERN), subido por el usuario 2026-10; facet_physics.csv",
    "ADY-IPAC14": "Ady & Kersevan, IPAC'14 WEPME038 / IPAC'19 TUPMP037 (TPMC)",
    "KIUBE-V04": "Geometría cámara KIUBE declarada en F65/F66 (supuesto_declarado_pendiente_fuente)",
    "F66-TPMC": "Mapa propio validado analíticamente (fase66/mapa_vacio_f66.json)",
}


def cargar_mapa_molflow(csv_path: str, xml_path: str,
                        n_bins: int = 20, r_pared_mm: float = 548.0):
    """Binnea presión de facetas de tapas por radio. Devuelve (r_mm, p_mbar)."""
    root = ET.parse(xml_path).getroot()
    verts = [(float(v.get("x")), float(v.get("y")), float(v.get("z")))
             for v in root.find("Geometry/Vertices")]
    import csv as _csv
    pres = []
    with open(csv_path) as fh:
        rd = _csv.DictReader(fh)
        for row in rd:
            pres.append(float(row["Pressure [mbar]"]))
    radios, zs = [], []
    for f in root.find("Geometry/Facets"):
        idx = [int(i.get("vertex")) for i in f.find("Indices")]
        vs = [verts[i] for i in idx]
        cx = sum(v[0] for v in vs) / len(vs)
        cy = sum(v[1] for v in vs) / len(vs)
        cz = sum(v[2] for v in vs) / len(vs)
        radios.append(np.hypot(cx, cy) * 10.0)  # cm -> mm
        zs.append(cz * 10.0)
    radios = np.array(radios)
    zs = np.array(zs)
    pres = np.array(pres)
    tapas = np.abs(zs) > 60.0
    bordes = np.linspace(0, r_pared_mm, n_bins + 1)
    r_c, p_c = [], []
    for b in range(n_bins):
        m = tapas & (radios >= bordes[b]) & (radios < bordes[b + 1])
        if m.any():
            r_c.append(float((bordes[b] + bordes[b + 1]) / 2))
            p_c.append(float(pres[m].mean()))
    return r_c, p_c


KB = 1.380649e-23


def balance_global(csv_path: str, xml_path: str, temperatura_k: float = 293.15,
                   q_total_mbar_l_s: float = 1.0023e-3):
    """Fracción absorbida = Q_bombas / Q_inyectado (chequeo de convergencia MC).

    Imping.rate de facet_physics.csv es por m²; se multiplica por área de
    faceta y por el sticking (0.25) de los puertos.
    """
    root = ET.parse(xml_path).getroot()
    verts = [np.array([float(v.get("x")), float(v.get("y")), float(v.get("z"))])
             for v in root.find("Geometry/Vertices")]
    import csv as _csv
    pres, zeta = [], []
    with open(csv_path) as fh:
        for row in _csv.DictReader(fh):
            zeta.append(float(row["Imping.rate"]))
    absorbido = 0.0
    for i, f in enumerate(root.find("Geometry/Facets")):
        st = float(f.find("Sticking").get("constValue"))
        if st < 0.05:
            continue
        idx = [int(iv.get("vertex")) for iv in f.find("Indices")]
        vs = [verts[k] for k in idx]
        area_cm2 = sum(0.5 * np.linalg.norm(np.cross(vs[k] - vs[0],
                                                     vs[k + 1] - vs[0]))
                       for k in range(1, len(vs) - 1))
        absorbido += zeta[i] * area_cm2 * 1e-4 * st  # moléculas/s
    q_abs_mbar_l_s = absorbido * KB * temperatura_k * 10.0
    return {"q_absorbida_mbar_l_s": q_abs_mbar_l_s,
            "q_inyectada_mbar_l_s": q_total_mbar_l_s,
            "fraccion_absorbida": q_abs_mbar_l_s / q_total_mbar_l_s}


def ejecutar(csv_path: str, xml_path: str, fase66_json: str,
             desorciones: int):
    balance = balance_global(csv_path, xml_path)
    r_mm, p_mbar = cargar_mapa_molflow(csv_path, xml_path)
    mapa66 = json.loads(Path(fase66_json).read_text())
    r66 = mapa66["mapa_presion"]["r_mm"]
    p66 = mapa66["mapa_presion"]["p_mbar"]
    p66_interp = np.interp(r_mm, r66, p66)
    p_arr = np.array(p_mbar)
    eps = 1e-12
    ratio = p_arr / np.maximum(p66_interp, eps)
    err_rel = np.abs(p_arr - p66_interp) / np.maximum(p66_interp, eps)

    perfil_molflow = PerfilPresion.desde_mapa(r_mm, p_mbar)
    perfil_propio = PerfilPresion.desde_mapa(r66, p66)
    casos = {}
    for nombre, perfil in (("mapa_molflow_f67", perfil_molflow),
                           ("mapa_tpmc_f66", perfil_propio)):
        res = simular_perdida_gas(
            corriente_ua=200.0, energia_final_mev=17.4,
            delta_e_kev_vuelta=80.0, b_media_t=1.77, perfil=perfil)
        casos[nombre] = {
            "fraccion_perdida": res.fraccion_perdida,
            "potencia_perdida_w": res.potencia_perdida_w,
            "estado": res.estado,
        }
    delta = (casos["mapa_molflow_f67"]["fraccion_perdida"]
             - casos["mapa_tpmc_f66"]["fraccion_perdida"])
    rel = (abs(delta) / casos["mapa_tpmc_f66"]["fraccion_perdida"]
           if casos["mapa_tpmc_f66"]["fraccion_perdida"] else None)

    salida = {
        "fase": 67,
        "version": "2.3.0",
        "fecha_utc": datetime.now(timezone.utc).isoformat(),
        "estado": "simulacion_teorica",
        "pq": "no_declarada_nunca",
        "desorciones_molflow": desorciones,
        "balance_global": balance,
        "mapa_molflow_f67": {"r_mm": r_mm, "p_mbar": p_mbar},
        "mapa_tpmc_f66_ref": {"r_mm": r66, "p_mbar": p66},
        "comparacion": {
            "ratio_medio_molflow_vs_propio": float(np.mean(ratio)),
            "error_rel_medio": float(np.mean(err_rel)),
            "error_rel_max": float(np.max(err_rel)),
            "perdida_stripping": casos,
            "delta_perdida_relativa": rel,
        },
        "fuentes": FUENTES,
        "pendientes": [
            "Convergencia Monte Carlo declarada (fracción absorbida registrada)",
            "Geometría KIUBE-V04 sigue siendo supuesto_declarado_pendiente_fuente",
            "Validación contra vacío medido: BLOQUEADA hasta haz real (memoria 3)",
        ],
        "nota_18_10": "Resultado informativo; no constituye decisión operacional.",
    }
    (EVID / "molflow_xval_f67.json").write_text(
        json.dumps(salida, indent=2, ensure_ascii=False))
    with open(EVID / "registro_vacio.jsonl", "a") as fh:  # append-only ALCOA+
        fh.write(json.dumps({
            "ts": salida["fecha_utc"], "fase": 67, "evento": "xval_molflow",
            "ratio_medio": salida["comparacion"]["ratio_medio_molflow_vs_propio"],
            "estado": "simulacion_teorica"}, ensure_ascii=False) + "\n")
    return salida


if __name__ == "__main__":
    import sys as _s
    out = ejecutar(_s.argv[1], _s.argv[2], _s.argv[3],
                   int(float(_s.argv[4])) if len(_s.argv) > 4 else 0)
    print(json.dumps(out["comparacion"], indent=2))
