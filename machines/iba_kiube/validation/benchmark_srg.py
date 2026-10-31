"""Benchmark de blindaje: Monte Carlo propio (Geant4, mc_benchmark/srg_bench) frente al SRG de IBA (MID 55672 Rev. G).
  Tabla 8: dosis a 1 m del ciclotrón dentro del búnker (0°, 90°, 180° respecto del blanco 18F), KIUBE 180 doble haz 2×100 µA.
  Tabla 9: dosis de contacto tras 200 / 220 / 230 cm de hormigón (2.35 g/cm³) frente al blanco.
Lee out/mc_srg/res/t<cm>_s<semilla>.txt y escribe validation/benchmark_srg.{json,png}.
Normalización de la fuente: BON21 (4.93e-3 n/p, FLUKA) y SRG (6.9e-3 n/p, IBA)."""
import glob, json, math, pathlib, re, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube import load_params                                                    # noqa: E402

P, _ = load_params()
RES = ROOT / "out/mc_srg/res"
geo = json.load(open(ROOT / "out/mc_srg/geometria.json"))
I_target_uA = 100.0                                   # SRG Tabla 8/9: KIUBE 180 doble haz «2 × 100 µA»
pps = I_target_uA * 1e-6 / 1.602176634e-19           # protones/s por blanco
Y = {"BON21": geo["yield_BON21"], "SRG": geo["yield_SRG"]}
S = {k: 2 * pps * y for k, y in Y.items()}           # n/s totales (2 blancos)
to_uSv_h = {k: s * 3600 * 1e-6 for k, s in S.items()}   # pSv/n → µSv/h

SRG8 = {"frente_0deg_interior": 1.2e8, "lateral_90deg_interior": 2.7e7, "opuesto_180deg_interior": 1.2e8}
SRG9 = {200: 3.62, 220: 0.79, 230: 0.40}


def leer(f):
    hdr = open(f).readline()
    N = int(re.search(r"N (\d+)", hdr).group(1))
    out = {}
    for l in open(f):
        if l.startswith("#"):
            continue
        n, mN, s2N, mG, s2G, *_ = l.split()
        out[n] = dict(N=N, mN=float(mN), s2N=float(s2N), mG=float(mG), s2G=float(s2G))
    return out


def combinar(lst, comp):
    Nt = sum(d["N"] for d in lst)
    k = "m" + comp; k2 = "s2" + comp
    s = sum(d[k] * d["N"] for d in lst); s2 = sum(d[k2] for d in lst)
    m = s / Nt
    var = max(s2 / Nt - m * m, 0.0) / Nt
    return m, math.sqrt(var)


runs = {}
for f in sorted(glob.glob(str(RES / "t*_s*.txt"))):
    t = int(re.search(r"t(\d+)_s", f).group(1))
    runs.setdefault(t, []).append(leer(f))
if not runs:
    sys.exit("sin resultados en out/mc_srg/res")
dets = list(next(iter(runs.values()))[0].keys())
res = {"fuente": dict(I_por_blanco_uA=I_target_uA, blancos=2, n_por_s=S, n_por_p=Y,
                      espectro="BON21 (FLUKA), isótropo", fisica="Geant4 11.4.3 QGSP_BIC_HP, G4NDL 4.7.1",
                      dosis="H*(10), ICRP 74", vr="importancia geométrica, capas de 10 cm, f = 3 por capa"),
       "espesores": {}, "tabla8": {}, "tabla9": {}}
for t, lst in sorted(runs.items()):
    res["espesores"][t] = {}
    for d in dets:
        sub = [x[d] for x in lst]
        mN, eN = combinar(sub, "N"); mG, eG = combinar(sub, "G")
        res["espesores"][t][d] = dict(n_pSv=mN, n_err=eN, g_pSv=mG, g_err=eG, tot_pSv=mN + mG,
                                      tot_err=math.hypot(eN, eG), N=sum(x["N"] for x in sub))
# Tabla 8: detectores interiores — promedio ponderado por varianza sobre todos los espesores (casi no dependen de t)
for d in list(SRG8) + ["a_1m_del_blanco_0deg"]:
    vals = [(r[d]["tot_pSv"], r[d]["tot_err"]) for r in res["espesores"].values() if r[d]["tot_err"] > 0]
    w = [1 / e ** 2 for _, e in vals]
    m = sum(v * wi for (v, _), wi in zip(vals, w)) / sum(w); e = 1 / math.sqrt(sum(w))
    res["tabla8"][d] = {"modelo_pSv_por_n": m, "err": e,
                        **{f"modelo_uSv_h_{k}": m * c for k, c in to_uSv_h.items()},
                        "SRG_uSv_h": SRG8.get(d), "cociente_modelo_SRG_conYsrg": (m * to_uSv_h["SRG"] / SRG8[d]) if d in SRG8 else None}
# Tabla 9: detector exterior central (1 m × 1 m, 20 cm de aire, contacto con el muro)
for t, r in res["espesores"].items():
    c = r["exterior_tras_muro"]
    res["tabla9"][t] = {"modelo_pSv_por_n": c["tot_pSv"], "err_rel": c["tot_err"] / c["tot_pSv"] if c["tot_pSv"] else None,
                        "fraccion_gamma": c["g_pSv"] / c["tot_pSv"] if c["tot_pSv"] else None,
                        **{f"modelo_uSv_h_{k}": c["tot_pSv"] * v for k, v in to_uSv_h.items()},
                        "SRG_uSv_h": SRG9.get(t),
                        "cociente_modelo_SRG_conYsrg": (c["tot_pSv"] * to_uSv_h["SRG"] / SRG9[t]) if t in SRG9 else None}
# décimo-espesor (TVL) del modelo y del SRG
# v0.27: ajuste ponderado (1/σ relativa) sobre TODOS los espesores corridos; con solo 200/220/230 un evento gamma
# de peso alto (220 cm, ±98 %) movía el TVL de 25 a 34 cm
ts = np.array(sorted(res["espesores"])); dX = "exterior_tras_muro"
hs = np.array([res["espesores"][t][dX]["tot_pSv"] for t in ts]); es = np.array([res["espesores"][t][dX]["tot_err"] for t in ts])
ok = (hs > 0) & (es < hs)                                   # descarta casos con error relativo ≥ 100 %
if ok.sum() >= 2:
    # corrección geométrica 1/r² desde el blanco (r = distancia blanco–cara exterior)
    r_out = np.array([P["sitio"]["vault_min"][1] / 2 + 100 + t * 10 + 100 - geo["sources"][0]["pos_mm"][1] for t in ts]) / 1000
    a, b = np.polyfit(ts[ok], np.log10(hs[ok] * r_out[ok] ** 2), 1, w=hs[ok] / es[ok])
    res["TVL_modelo_cm"] = -1 / a
    res["TVL_ajuste"] = dict(espesores_cm=[int(t) for t in ts[ok]], ponderado="1/σ_rel", nota="corregido por 1/r²")
tS = np.array(sorted(SRG9)); hS = np.array([SRG9[t] for t in tS])
r_S = np.array([2.0 + 0.1 + t / 100 + 0.1 - geo["sources"][0]["pos_mm"][1] / 1000 for t in tS])
aS, _ = np.polyfit(tS, np.log10(hS * r_S ** 2), 1)
res["TVL_SRG_cm"] = -1 / aS
json.dump(res, open(ROOT / "validation/benchmark_srg.json", "w"), indent=1, ensure_ascii=False, default=float)

# ---------------------------------------------------------------- figura
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, mk, col in (("SRG", "o", "#c0392b"), ("BON21", "s", "#2c6fbb")):
    y = np.array([res["tabla9"][t][f"modelo_uSv_h_{k}"] for t in ts])
    e = np.array([res["tabla9"][t]["err_rel"] or 0 for t in ts]) * y
    ax[0].errorbar(ts, y, yerr=e, fmt=mk + "-", color=col, capsize=3, label=f"Modelo Geant4 (fuente {k}: {Y[k]:.2e} n/p)")
ax[0].plot(tS, hS, "k^", ms=10, label="SRG Tabla 9 (IBA)")
ax[0].set_yscale("log"); ax[0].set_xlabel("Espesor del muro de hormigón 2.35 g/cm³ [cm]")
ax[0].set_ylabel("H*(10) de contacto frente al blanco [µSv/h]"); ax[0].grid(alpha=.3, which="both")
ax[0].set_title(f"Tras el muro — KIUBE 180, 2 × 100 µA\nTVL modelo {res.get('TVL_modelo_cm', float('nan')):.1f} cm · SRG {res['TVL_SRG_cm']:.1f} cm")
ax[0].legend(fontsize=8)
lab = {"frente_0deg_interior": "0° (frente)", "lateral_90deg_interior": "90°", "opuesto_180deg_interior": "180°"}
xs = np.arange(3)
mm = [res["tabla8"][d]["modelo_uSv_h_SRG"] for d in lab]; ee = [res["tabla8"][d]["err"] * to_uSv_h["SRG"] for d in lab]
ax[1].bar(xs - 0.2, mm, 0.4, yerr=ee, color="#c0392b", capsize=3, label="Modelo (fuente SRG), junto al muro")
ax[1].bar(xs + 0.2, [SRG8[d] for d in lab], 0.4, color="#555", label="SRG Tabla 8 «a 1 m»")
ax[1].set_xticks(xs); ax[1].set_xticklabels(lab.values()); ax[1].set_yscale("log")
ax[1].set_ylabel("H*(10) dentro del búnker [µSv/h]"); ax[1].set_title("Dentro del búnker (Tabla 8)"); ax[1].legend(fontsize=8)
ax[1].grid(alpha=.3, which="both", axis="y")
fig.tight_layout(); fig.savefig(ROOT / "validation/benchmark_srg.png", dpi=150)
print(json.dumps({"tabla8": res["tabla8"], "tabla9": res["tabla9"], "TVL_modelo_cm": res.get("TVL_modelo_cm"),
                  "TVL_SRG_cm": res["TVL_SRG_cm"]}, indent=1, ensure_ascii=False, default=float))
