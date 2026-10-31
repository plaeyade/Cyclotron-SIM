"""Barrido del gap fuente–puller de la región central del KIUBE (preset 'kiube', v0.23).

Variables: gap_sp (gap chimenea–puller), radio de la fuente y regla de azimut de la fuente:
  'modelo' : la chimenea queda (chimenea + gap + puller)/r antes del borde de entrada del dee (puller pegado al borde)
  'IBA40'  : la chimenea se fija 40° antes del borde de entrada, como en el mapa Opera de IBA (TUD03T diap. 5);
             src_dtheta compensa el cambio de gap
Métricas (mallas gruesas, 3 vueltas, como optimize_sources): descentrado de la órbita, fases que sobreviven, ganancia
de energía por vuelta y campo en el gap fuente–puller (V_dee / gap).
Uso: python fem/central_region_barrido_gap.py  → out/central_region/barrido_gap/resultado.json + figura
"""
from __future__ import annotations
import itertools, json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import central_region as cr
from kiube import load_params

GAPS = [1.5, 2.0, 2.6]
RADIOS = [12.5, 14.0]
REGLAS = ["modelo", "IBA40"]


def caso(pr0, P, gap, r, regla, fases=(-60.0, 20.0, 10.0), vueltas=3.0):
    q = dict(pr0, h_gap=0.35, h_nose=1.6, h_mid=8.0, h_far=30.0, R_dom=min(pr0["R_dom"], 220.0), R_fino=40.0)
    q["gap_sp"] = gap
    reach = q["chimney_r"] + gap + q["puller_len"]
    dth = 0.0 if regla == "modelo" else math.degrees(reach / r) - 40.0
    q["fuentes"] = [dict(F, src_r=r, src_dtheta=dth) for F in pr0["fuentes"]]
    mesh, gfu, info, st = cr.solve(q, None)
    Ef = cr.FEField(mesh, gfu, q.get("mitad", False))
    F = info["fuentes"][0]
    res, *_ = cr.run_species(Ef, q, F, P, fases, vueltas, 0, R_stop=0.9 * q["R_dom"], beam=False)
    ok = [x for x in res["barrido"] if x["sobrevive"] and x["centro_mm"]]
    dc = min((math.hypot(x["centro_mm"][0], x["centro_mm"][1]) for x in ok), default=None)
    Tm = max((x["T_final_keV"] for x in ok), default=None)
    return dict(gap_sp=gap, src_r=r, regla=regla, src_dtheta=round(dth, 2), angulo_antes_borde=round(math.degrees(reach / r) - dth, 1),
                fases_ok=len(ok), fases_total=len(res["barrido"]), desc_centro_mm=dc, T_max_keV_3v=Tm,
                E_gap_MV_m=q["V_dee"] / gap / 1e3, malla=st["ne"])


if __name__ == "__main__":
    P, _ = load_params()
    pr0 = dict(cr.PRESETS["kiube"])
    out = ROOT / "out/central_region/barrido_gap"; out.mkdir(parents=True, exist_ok=True)
    filas = []
    for gap, r, regla in itertools.product(GAPS, RADIOS, REGLAS):
        f = caso(pr0, P, gap, r, regla)
        print("[gap]", json.dumps(f), flush=True)
        filas.append(f)
        json.dump(filas, open(out / "resultado.json", "w"), indent=1, default=float)
