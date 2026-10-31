"""Campo eléctrico máximo en la región fuente–puller (v0.31): |E| en una rejilla fina delante de la chimenea y de las
mordazas del puller, a V_dee = 40 kV, frente al criterio de Kilpatrick (8.3 MV/m a 40.6 MHz; fem/README_RF.md).
Uso: python fem/central_region_campo_max.py tag [clave=valor ...]  → out/central_region/campo_max/<tag>.json
"""
import json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))
import central_region as cr
KP = 8.3                                                       # MV/m (Kilpatrick a 40.6 MHz)

if __name__ == "__main__":
    tag = sys.argv[1]
    pr = dict(cr.PRESETS["kiube"])
    for kv in sys.argv[2:]:
        k, v = kv.split("="); pr[k] = type(pr[k])(v)
    mesh, gfu, info, st = cr.solve(pr, None)
    Ef = cr.FEField(mesh, gfu, False)
    F = info["fuentes"][0]; S = np.array(F["S"]); en = np.array(F["e_n"]); et = np.array([-en[1], en[0]])
    rc, gap = pr["chimney_r"], pr["gap_sp"]
    s_ = np.arange(rc - 0.1, rc + gap + 1.5, 0.05); u_ = np.arange(-pr["puller_w"] / 2 - 0.5, pr["puller_w"] / 2 + 0.51, 0.1)
    z_ = np.arange(-8.0, 8.01, 0.05)
    best = []
    for u in u_:
        SS, ZZ = np.meshgrid(s_, z_, indexing="ij")
        P = np.c_[(S[None, :] + SS.ravel()[:, None] * en[None, :] + u * et[None, :]), ZZ.ravel()]
        E = Ef(P); ok = Ef._last_ok
        m = np.linalg.norm(E, axis=1) * pr["V_dee"] * 1e3 / 1e6        # V/mm por voltio → MV/m
        m[~ok] = 0
        i = int(np.argmax(m)); best.append((m[i], float(SS.ravel()[i]), float(u), float(ZZ.ravel()[i])))
    mx = max(best)
    # campo en el eje del gap (u = z = 0), punto medio
    Pm = np.array([[*(S + (rc + gap / 2) * en), 0.0]])
    Em = float(np.linalg.norm(Ef(Pm)) * pr["V_dee"] * 1e3 / 1e6)
    out = dict(tag=tag, sets=sys.argv[2:], E_max_MV_m=mx[0], donde_s_u_z_mm=mx[1:], E_max_sobre_Kp=mx[0] / KP,
               E_eje_gap_MV_m=Em, malla=st["ne"], nota="rejilla de 0.05 mm en vacío; el máximo real en una arista viva es singular (depende de la malla)")
    d = ROOT / "out/central_region/campo_max"; d.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(d / f"{tag}.json", "w"), indent=1, default=float)
    print(json.dumps(out, default=float))
