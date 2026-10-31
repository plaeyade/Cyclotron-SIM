"""Región central v0.31: forma de las mordazas del puller (bisel del borde interior y escalón de la abertura).
Entradas: out/central_region/puller/v31_*, v29b_base; out/central_region/menisco*/transmision_jp1440*.json;
out/central_region/campo_max/*.json.  Salidas: validation/region_central_v31.{png,json}
Métrica comparable entre casos: fracción de la corriente continua de la fuente que llega a 6 vueltas
= transmisión × (ancho de la ventana de lanzamiento)/360°, porque cada caso lanza en su propia ventana aceptada ±10°."""
import json, math, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; D = ROOT / "out/central_region/puller"

def fila(tag):
    e = json.load(open(D / tag / "resultado.json"))["especies"][0]
    n = sum(e["perdidas_haz"].values()); p = e["perdidas_haz"]["vive"] / n
    lo, hi = e["ventana_fase_aceptada"]; win = hi - lo + 20.0
    return dict(tag=tag, T=100 * p, err=100 * math.sqrt(p * (1 - p) / n), ventana_deg=win, dc=100 * p * win / 360,
                dc_err=100 * math.sqrt(p * (1 - p) / n) * win / 360, z=100 * e["perdidas_haz"]["z"] / n, metal=100 * e["perdidas_haz"]["metal"] / n)

ch = [(0, "v29b_base"), (1, "v31_ch1"), (2, "v31_ch2"), (3, "v31_ch3"), (4, "v31_ch4"), (5, "v31_ch5"), (6, "v31_ch6")]
F = [dict(c=c, **fila(t)) for c, t in ch]
otros = {t: fila(t) for t in ("v31_ch2a30", "v31_ch2a60", "v31_ch3a30", "v31_e7", "v31_e8", "v31_l8e12", "v31_p8ch2", "v31_s12ch2")
         if (D / t / "resultado.json").exists()}
men = {}
for lab, f in (("sin bisel, menisco", "menisco/transmision_jp1440.json"), ("sin bisel, plasma plano", "menisco/transmision_jp1440_plano.json"),
               ("bisel 5 mm, menisco", "menisco_ch5/transmision_jp1440.json"), ("bisel 5 mm, plasma plano", "menisco_ch5/transmision_jp1440_plano.json")):
    R = json.load(open(ROOT / "out/central_region" / f)); men[lab] = dict(T=R["transmision_pct"], err=R["error_pct"], dc=R["transmision_pct"] * 120 / 360,
                                                                         perdidas=R["perdidas_pct"])
E = {t: json.load(open(ROOT / f"out/central_region/campo_max/{t}.json")) for t in ("base", "ch5")}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
ax[0].errorbar([f["c"] for f in F], [f["dc"] for f in F], [f["dc_err"] for f in F], fmt="o-", color="#16a085", capsize=3, label="45°")
for t, x, mk in (("v31_ch2a30", 2, "s"), ("v31_ch3a30", 3, "s"), ("v31_ch2a60", 2, "^")):
    if t in otros:
        ax[0].errorbar([x + 0.08], [otros[t]["dc"]], [otros[t]["dc_err"]], fmt=mk, color="#8e44ad", capsize=3,
                       label=("30°" if "a30" in t else "60°") if t != "v31_ch3a30" else None)
ax[0].axvline(5, color="#d0182a", lw=1); ax[0].text(5.1, 9.2, "elegido", color="#d0182a", fontsize=8)
ax[0].set_xlabel("bisel del borde interior de las mordazas [mm]"); ax[0].set_ylabel("% de la corriente DC de la fuente a 6 vueltas")
ax[0].set_title("Bisel (plasma plano, 400 partículas)"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
labs = list(men); x = np.arange(len(labs))
ax[1].bar(x, [men[l]["T"] for l in labs], yerr=[men[l]["err"] for l in labs], capsize=4,
          color=["#16a085" if "menisco" in l else "#888" for l in labs])
for i, l in enumerate(labs):
    ax[1].text(i, 2, f"z {men[l]['perdidas']['z']:.0f} %\nmetal {men[l]['perdidas']['metal']:.0f} %", ha="center", fontsize=7, color="w")
ax[1].set_xticks(x); ax[1].set_xticklabels(labs, fontsize=7, rotation=12); ax[1].set_ylabel("transmisión 6 vueltas [%] (−80…+40°)")
ax[1].set_title("Confirmación con el menisco autoconsistente"); ax[1].grid(alpha=.3, axis="y")
ax[2].bar([0, 1], [E["base"]["E_max_MV_m"], E["ch5"]["E_max_MV_m"]], color=["#888", "#16a085"])
ax[2].axhline(8.3, color="#d0182a", ls="--", lw=1); ax[2].text(0.6, 8.8, "Kilpatrick 8.3 MV/m (40.6 MHz)", fontsize=7, color="#d0182a")
ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["sin bisel\n(arista viva)", "bisel 5 mm"]); ax[2].set_ylabel("|E| máx. en la región fuente–puller [MV/m]")
ax[2].set_title("Campo máximo a 40 kV"); ax[2].grid(alpha=.3, axis="y")
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_v31.png", dpi=130)
json.dump(dict(bisel_45=F, otros=otros, menisco=men, campo_max={k: dict(E_max_MV_m=v["E_max_MV_m"], E_max_sobre_Kp=v["E_max_sobre_Kp"],
                                                                     donde_s_u_z_mm=v["donde_s_u_z_mm"]) for k, v in E.items()},
               elegido=dict(puller_chaflan=5.0, puller_chaflan_ang=45.0)),
          open(ROOT / "validation/region_central_v31.json", "w"), indent=1, ensure_ascii=False, default=float)
print([(f["c"], round(f["dc"], 2)) for f in F])
