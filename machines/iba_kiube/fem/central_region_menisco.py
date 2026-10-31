"""Región central con el menisco autoconsistente del plasma (v0.30).

Para cada nivel de voltaje de fem/menisco.py (V = V_dee·cos φ en la emisión) se resuelve el campo global con el plasma
detrás de ese menisco (penalización) y se siguen las partículas emitidas en las fases cuyo cos φ cae en ese tramo. Las
partículas arrancan sobre el menisco. Mismas fases (−80…+40°), número (400) y criterio que el caso base.
Uso: python fem/central_region_menisco.py [--jp 1440]  → out/central_region/menisco/transmision_jp<jp>.json
"""
import argparse, json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import central_region as cr
from kiube import load_params

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jp", type=int, default=1440)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--plano", action="store_true", help="control: misma cavidad honda con el plasma plano en la cara original")
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = pathlib.Path(a.out) if a.out else ROOT / "out/central_region/menisco"
    M = json.load(open(out / f"menisco_jp{a.jp}.json"))
    niv = sorted(M["niveles"], key=lambda r: -r["cos"])
    cs = np.array([r["cos"] for r in niv])
    bordes = [(cs[i] + cs[i + 1]) / 2 for i in range(len(cs) - 1)]          # cos decreciente
    P, _ = load_params()
    pr0 = dict(cr.PRESETS["kiube"]); pr0["_nombre"] = "kiube"; pr0["puller_slit_h"] = 10.0
    pr0["menisco_s_cav"] = pr0["chimney_r"] - 1.2
    for kv in a.set:
        k, v = kv.split("="); pr0[k] = type(pr0[k])(v)
    rng = np.random.default_rng(7)
    PH = rng.uniform(-80.0, 40.0, a.n)                                       # = ventana del caso base (−70…30 ± 10)
    cphi = np.cos(np.radians(PH))
    k_of = np.array([int(np.sum(c < np.array(bordes))) for c in cphi])      # índice de nivel por partícula
    filas, tot = [], dict(vive=0, z=0, metal=0, r=0, nan=0)
    for k, r in enumerate(niv):
        sel = np.flatnonzero(k_of == k)
        if not len(sel):
            continue
        pr = dict(pr0, _menisco=([pr0["plasma_depth"], 0, 0, 0, 0] if a.plano else r["coef"]), _PH=PH[sel].tolist())
        mesh, gfu, info, st = cr.solve(pr, None)
        Ef = cr.FEField(mesh, gfu, False)
        F = info["fuentes"][0]
        res, trs, B, Bf = cr.run_species(Ef, pr, F, P, (-70.0, 30.0, 5.0), 6.0, len(sel), 0.9 * pr["R_dom"])
        ph = res["perdidas_haz"]
        for m_, v in ph.items():
            tot[m_] += v
        filas.append(dict(cos=r["cos"], n=int(len(sel)), perdidas=ph, coef=r["coef"], aceptancia_deg=res["ancho_aceptancia_deg"]))
        print(f"[menisco-global] cos {r['cos']}: {len(sel)} partículas → {ph}", flush=True)
        del mesh, gfu, Ef
    n = sum(tot.values()); p = (tot["vive"] + tot["r"]) / n
    R = dict(jp_A_m2=a.jp, plano=a.plano, n=n, transmision_pct=100 * p, error_pct=100 * math.sqrt(p * (1 - p) / n),
             perdidas_pct={k: 100 * v / n for k, v in tot.items()}, por_nivel=filas)
    json.dump(R, open(out / f"transmision_jp{a.jp}{'_plano' if a.plano else ''}.json", "w"), indent=1, default=float)
    print("[menisco-global]", json.dumps({k: v for k, v in R.items() if k != "por_nivel"}, default=float))
