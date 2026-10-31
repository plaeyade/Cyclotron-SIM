"""Figuras de los barridos de la región central (v0.23 gap, v0.24 puller).
Entradas: validation/region_central_gap.json y out/central_region/puller/barrido.json
Salidas: validation/region_central_gap.png, validation/region_central_puller.png y validation/region_central_puller.json"""
import json, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
TECHO = 105 / 120 * 100            # las partículas se lanzan en la ventana aceptada ±10° (≈ 120°); aceptancia ≈ 105°

# ---------- gap fuente–puller (v0.23)
g = json.load(open(ROOT / "validation/region_central_gap.json"))
G, T, S = np.array(g["gap_mm"]), np.array(g["transmision_pct"]), np.array(g["error_pct"])
c = np.polyfit(G, T, 2, w=1 / S)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.errorbar(G, T, S, fmt="o", color="#16a085", capsize=4, label="transmisión (400 partículas, 6 vueltas, puller 5 mm)")
gg = np.linspace(1.9, 3.3, 100); ax.plot(gg, np.polyval(c, gg), "--", color="#16a085", lw=1, label=f"ajuste: óptimo {g['optimo_gap_mm']:.2f} mm")
ax.axhline(TECHO, color="#888", ls=":", lw=1, label="techo por aceptancia de fase (≈ 105°/120°)")
ax.axvline(2.8, color="#8e44ad", lw=1); ax.text(2.82, 5, "v0.23: 2.8 mm", fontsize=7, color="#8e44ad")
ax.set_xlabel("gap fuente–puller [mm]"); ax.set_ylabel("%"); ax.set_ylim(0, 95); ax.grid(alpha=.3); ax.legend(fontsize=7.5)
ax.set_title("Gap fuente–puller (v0.23)")
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_gap.png", dpi=140)

# ---------- puller (v0.24)
P = json.load(open(ROOT / "out/central_region/puller/barrido.json"))
ref = json.load(open(ROOT / "out/central_region/kiube/resultado.json"))["especies"][0] if False else None
sl = sorted([f for f in P if set(f["sets"]) == {"puller_slit_h"}], key=lambda f: float(f["sets"]["puller_slit_h"]))
# referencia v0.23 (abertura 5 mm): 19.5 ± 2.0 %, metal 52.5 %, z 28 %
H = [5.0] + [float(f["sets"]["puller_slit_h"]) for f in sl]
TT = [19.5] + [f["transmision_pct"] for f in sl]; EE = [2.0] + [f["error_pct"] for f in sl]
MM = [52.5] + [f["perdidas"]["metal"] for f in sl]; ZZ = [28.0] + [f["perdidas"]["z"] for f in sl]
o = np.argsort(H); H, TT, EE, MM, ZZ = (np.array(v)[o] for v in (H, TT, EE, MM, ZZ))
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
ax[0].errorbar(H, TT, EE, fmt="o-", color="#16a085", capsize=4, label="transmisión del haz")
ax[0].axhline(TECHO, color="#888", ls=":", lw=1, label="techo por aceptancia de fase")
ax[0].axvline(10, color="#8e44ad", lw=1); ax[0].text(10.1, 5, "v0.24: 10 mm", fontsize=7, color="#8e44ad")
ax[0].axvline(6, color="#888", lw=0.8, ls="--"); ax[0].text(6.1, 80, "rendija de\nla chimenea", fontsize=7, color="#666")
ax[0].set_xlabel("abertura vertical del puller [mm]"); ax[0].set_ylabel("%"); ax[0].set_ylim(0, 95); ax[0].grid(alpha=.3)
ax[0].legend(fontsize=7.5); ax[0].set_title("Puller: abertura vertical (gap 2.8, ancho 7, prof. 3)")
ax[1].plot(H, MM, "s-", color="#d0182a", label="choca con metal (mordazas del puller, chimenea)")
ax[1].plot(H, ZZ, "^-", color="#3b6fd8", label="se pierde en vertical (|z| > 7.8 mm)")
ax[1].set_xlabel("abertura vertical del puller [mm]"); ax[1].set_ylabel("% del haz"); ax[1].grid(alpha=.3)
ax[1].legend(fontsize=7.5); ax[1].set_title("Dónde se pierde el haz")
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_puller.png", dpi=140)
otros = [dict(tag=f["tag"], sets=f["sets"], transmision_pct=f["transmision_pct"], error_pct=f["error_pct"],
              aceptancia_deg=f["aceptancia_deg"], perdidas=f["perdidas"]) for f in P if set(f["sets"]) != {"puller_slit_h"}]
json.dump(dict(abertura_mm=H.tolist(), transmision_pct=TT.tolist(), error_pct=EE.tolist(), metal_pct=MM.tolist(),
               z_pct=ZZ.tolist(), otros=otros, elegido=dict(puller_slit_h=10, puller_w=7, puller_len=3, gap_sp=2.8, ap_half=8),
               techo_pct=TECHO, nota="preset kiube, malla fina, 6 vueltas, 400 partículas, fuente a r 12.5 y 40.3° del borde"),
          open(ROOT / "validation/region_central_puller.json", "w"), indent=1, default=float)
print("ok")

# ---------- rendija de la chimenea (v0.25): transmisión y corriente útil ∝ altura × transmisión
R = [f for f in P if f["tag"].startswith("rend")] + [f for f in P if f["tag"] == "slit10"]
hh = np.array([float(f["sets"].get("slit_h", 6)) for f in R]); o = np.argsort(hh); R = [R[i] for i in o]; hh = hh[o]
tr = np.array([f["transmision_pct"] for f in R]); er = np.array([f["error_pct"] for f in R])
cu = hh * tr / 100; ce = hh * er / 100
fig, ax = plt.subplots(figsize=(6.6, 4.4))
ax.errorbar(hh, tr, er, fmt="o-", color="#16a085", capsize=4, label="transmisión [%]")
a2 = ax.twinx(); a2.errorbar(hh, cu, ce, fmt="s--", color="#d0182a", capsize=4, label="corriente útil ∝ altura × transmisión")
a2.set_ylabel("corriente útil [u. rel.]", color="#d0182a"); a2.set_ylim(0, 1.3 * cu.max())
ax.axvline(6, color="#8e44ad", lw=1); ax.text(6.1, 5, "v0.25: 6 mm", fontsize=7, color="#8e44ad")
ax.set_xlabel("altura de la rendija de la chimenea [mm]"); ax.set_ylabel("transmisión [%]", color="#16a085"); ax.set_ylim(0, 70)
ax.grid(alpha=.3); ax.set_title("Rendija de la fuente PIG (puller 10 × 7 × 3, gap 2.8)")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper right")
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_rendija.png", dpi=140)
json.dump(dict(slit_h_mm=hh.tolist(), transmision_pct=tr.tolist(), error_pct=er.tolist(), corriente_rel=cu.tolist(),
               elegido_mm=6.0, criterio="máxima corriente útil (la emisión de la PIG ∝ área de la rendija)"),
          open(ROOT / "validation/region_central_rendija.json", "w"), indent=1, default=float)
print("rendija ok")
