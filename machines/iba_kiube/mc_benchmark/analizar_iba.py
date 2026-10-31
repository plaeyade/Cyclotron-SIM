"""Análisis del MC del búnker IBA (RPOS F89). Para cada cara usa SOLO las corridas sesgadas hacia esa cara:
perfil de H*(10) medio sobre la cara (proyección de la bóveda) y máximo en parches de 50×50 cm por capa de 10 cm;
ajuste exponencial; exterior directo con los espesores IBA; espesor requerido por objetivo (300 µA, 2000 h/año).
Salida: out/mc_iba/resultado_iba.json y validation/iba_mapa.png"""
import glob, json, math, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from leer_iba import leer
ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "out/mc_iba/res"
# caras en coordenadas mundo: (eje, signo, cara interior mm, espesor IBA mm, rango perpendicular a promediar)
FL, CE = -1120.0, 1380.0
CARAS = {
    "este": dict(ax=1, sg=1, p0=2000, t=2200, perp={"x": (-1970, 1970), "z": (FL, CE)}),
    "oeste": dict(ax=1, sg=-1, p0=2000, t=2400, perp={"x": (-1970, 1970), "z": (FL, CE)}),
    "sur": dict(ax=0, sg=1, p0=1970, t=2500, perp={"y": (-2000, 2000), "z": (FL, CE)}),
    "techo": dict(ax=2, sg=1, p0=1380, t=2300, perp={"x": (-1970, 1970), "y": (-2000, 2000)}),
    "norte": dict(ax=0, sg=-1, p0=1970, t=7240, perp={"y": (-4400, 4200), "z": (FL, CE)}),
}
# objetivos (µSv/h con haz) — F88: restricción 0.3 mSv/año público, 2 mSv/año personal expuesto; 2000 h/año
OBJ = {"este": (1.0, "radiofarmacia, personal expuesto, T=1, 2 mSv/año"), "oeste": (0.75, "exterior/circulación T=1/5, 0.3 mSv/año"),
       "sur": (6.0, "área verde T=1/40, 0.3 mSv/año"), "techo": (6.0, "área verde vacía sobre el techo T=1/40, 0.3 mSv/año"),
       "norte": (0.15, "control/enfriamiento tras el laberinto T=1, 0.3 mSv/año")}


def cargar(cara):
    fs = sorted(glob.glob(str(D / f"{cara}_s*_malla.bin")))
    Ms = [leer(f) for f in fs]
    return Ms


def perfil(Ms, c):
    M0 = Ms[0]; coords = {"x": M0["x"], "y": M0["y"], "z": M0["z"]}; ax = "xyz"[c["ax"]]
    sel = [np.ones_like(coords[k], bool) for k in "xyz"]
    for k, (a, b) in c["perp"].items():
        i = "xyz".index(k); sel[i] = (coords[k] > a) & (coords[k] < b)
    dcoord = c["sg"] * coords[ax] - c["p0"]                     # profundidad desde la cara interior [mm]
    res = []
    for kk in range(-2, int((c["t"] + 1000) / 100)):
        m = (dcoord > kk * 100) & (dcoord < (kk + 1) * 100)
        if not m.any():
            continue
        s = list(sel); s[c["ax"]] = s[c["ax"]] & m
        vals, maxs = [], []
        for M in Ms:
            H = (M["Hn"] + M["Hg"])[np.ix_(s[0], s[1], s[2])]
            vals.append(H.mean())
            # parches de 5×5 vóxeles (50×50 cm) en el plano perpendicular
            Hp = H.squeeze()
            if Hp.ndim == 2 and min(Hp.shape) >= 5:
                a, b = Hp.shape[0] // 5 * 5, Hp.shape[1] // 5 * 5
                P = Hp[:a, :b].reshape(a // 5, 5, b // 5, 5).mean(axis=(1, 3)); maxs.append(P)
        v = np.array(vals); n = len(v)
        pm = np.mean(maxs, axis=0).max() if maxs else float("nan")
        res.append(dict(prof_cm=(kk + 0.5) * 10, H=float(v.mean()), se=float(v.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan"), parche_max=float(pm)))
    return res


def ajustar(P, t_cm, lo=40, hi=None):
    hi = hi or t_cm - 10
    z = np.array([p["prof_cm"] for p in P]); H = np.array([p["H"] for p in P]); se = np.array([p["se"] for p in P])
    ok = (z >= lo) & (z <= hi) & (H > 0) & np.isfinite(se) & (se < 0.5 * H)
    if ok.sum() < 4:
        return None
    w = 1 / np.maximum(se[ok] / H[ok], 0.05)
    b, a = np.polyfit(z[ok], np.log10(H[ok]), 1, w=w)
    return dict(a=float(a), b=float(b), TVL_cm=float(-1 / b), n_puntos=int(ok.sum()), rango_cm=[lo, hi])


def main():
    out = {"corridas": {}, "caras": {}}
    for cara, c in CARAS.items():
        Ms = cargar(cara)
        if not Ms:
            continue
        out["corridas"][cara] = dict(lotes=len(Ms), historias=float(sum(M["N"] for M in Ms)))
        P = perfil(Ms, c); t = c["t"] / 10
        fit = ajustar(P, t)
        ext = [p for p in P if t < p["prof_cm"] < t + 30]
        H_ext = float(np.mean([p["H"] for p in ext])) if ext else None
        se_ext = float(np.sqrt(np.sum([p["se"] ** 2 for p in ext])) / len(ext)) if ext else None
        pm_ext = float(np.max([p["parche_max"] for p in ext])) if ext else None
        obj, txt = OBJ[cara]
        r = dict(espesor_IBA_cm=t, perfil=P, ajuste=fit, exterior_0_30cm_uSv_h=H_ext, exterior_se=se_ext, exterior_parche_max_uSv_h=pm_ext,
                 objetivo_uSv_h=obj, objetivo=txt)
        if fit:
            # espesor requerido usando la cara media + factor parche/medio medido en la última capa interior con estadística
            ratio = [p["parche_max"] / p["H"] for p in P if fit["rango_cm"][0] <= p["prof_cm"] <= fit["rango_cm"][1] and p["H"] > 0 and np.isfinite(p["parche_max"])]
            fp = float(np.median(ratio)) if ratio else 1.0
            treq = (math.log10(obj / fp) - fit["a"]) / fit["b"]
            r.update(factor_parche=fp, espesor_requerido_cm=round(treq, 0), modulo_13cm=13 * math.ceil(treq / 13),
                     ahorro_vs_IBA_cm=round(t - 13 * math.ceil(treq / 13), 0),
                     exterior_ajuste_uSv_h=float(10 ** (fit["a"] + fit["b"] * (t + 5))))
        # método directo (primario): menor profundidad z tal que el máximo de parche 50×50 cm es ≤ objetivo en z y en
        # todas las capas más profundas dentro del muro; interpolación log-lineal entre capas. El ajuste global queda como chequeo.
        Pin = [p for p in P if 0 < p["prof_cm"] < t and np.isfinite(p["parche_max"]) and p["parche_max"] > 0]
        zr = None
        for i, p in enumerate(Pin):
            if all(q["parche_max"] <= obj for q in Pin[i:]):
                if i == 0:
                    zr = p["prof_cm"]
                else:
                    a, b = Pin[i - 1], p
                    la, lb = math.log10(a["parche_max"]), math.log10(b["parche_max"])
                    zr = a["prof_cm"] + (math.log10(obj) - la) / (lb - la) * (b["prof_cm"] - a["prof_cm"])
                break
        if zr is not None:
            zr += 5                                                    # centro de capa → cara de salida de la capa
            r.update(espesor_requerido_directo_cm=round(zr, 0), modulo_13cm_directo=13 * math.ceil(zr / 13),
                     ahorro_directo_vs_IBA_cm=round(t - 13 * math.ceil(zr / 13), 0))
        out["caras"][cara] = r
    json.dump(out, open(ROOT / "out/mc_iba/resultado_iba.json", "w"), indent=1, default=float)
    for k, r in out["caras"].items():
        print(k, out["corridas"][k], "TVL", r["ajuste"] and round(r["ajuste"]["TVL_cm"], 1), "ext", r["exterior_0_30cm_uSv_h"], "±", r["exterior_se"],
              "parche", r["exterior_parche_max_uSv_h"], "req", r.get("espesor_requerido_cm"), "mod", r.get("modulo_13cm"), "ahorro", r.get("ahorro_vs_IBA_cm"),
              "| directo", r.get("espesor_requerido_directo_cm"), r.get("modulo_13cm_directo"), r.get("ahorro_directo_vs_IBA_cm"))
    return out


if __name__ == "__main__":
    main()
