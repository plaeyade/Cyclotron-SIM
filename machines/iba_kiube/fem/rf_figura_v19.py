"""Figura de v0.19: espectro de modos con y sin puente central y curva del sintonizador (out/fem_rf/resultado_*.json)."""
import json, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
R = lambda t: json.load(open(ROOT / f"out/fem_rf/resultado_{t}.json"))                       # noqa: E731
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
casos = [("v0.18: dees sin acoplar", "v18_z75", "#888888"), ("v0.19: puente central", "v19_solo_puente", "#16a085"),
         ("sin puente + pistón en el dee 1", "v19_sin_puente_con_sint", "#d0182a")]
for k, (lab, t, c) in enumerate(casos):
    try:
        for i, m in enumerate(sorted(R(t)["modos"], key=lambda m: m["f_MHz"])):
            f = m["f_MHz"]; ax[0].plot([f, f], [k + 0.6, k + 1.3], color=c, lw=3)
            nom = m["fase"].split(" (")[0]
            if nom == "otro":
                nom = "solo dee 0" if m["V1_sobre_V0"] < 0.1 else "solo dee 1" if m["V1_sobre_V0"] > 10 else "otro"
            if t == "v18_z75":
                nom = "en fase y contrafase" if i == 0 else ""
                if i: continue
            ax[0].text(f, k + 1.33 + 0.22 * (i % 2), f"{f:.1f} {nom}", ha="center", fontsize=7)
    except FileNotFoundError:
        pass
ax[0].set_yticks([1, 2, 3]); ax[0].set_yticklabels([c[0] for c in casos], fontsize=8)
ax[0].set_xlabel("frecuencia [MHz]"); ax[0].set_xscale("log"); ax[0].set_xlim(30, 300); ax[0].set_ylim(0.4, 3.9)
ax[0].set_title("Modos de dee: el puente aísla el modo acelerador"); ax[0].grid(alpha=.3, axis="x")
pts = []
for g in (9999, 200, 80):
    try:
        m = R(f"v19s_t40_g{g}")["modos"][0]
        depth = 0 if g == 9999 else 630 - (75 + g)
        pts.append((depth, m["f_MHz"], m["P_kW_40kV"]))
    except FileNotFoundError:
        pass
if pts:
    d, f, p = zip(*pts); f0 = f[0]
    ax[1].plot(d, [(x - f0) * 1e3 / 2 for x in f], "o-", ms=6, color="#3b6fd8", label="Δf por pistón [kHz] (mitad del par simétrico del FE)")
    ax2 = ax[1].twinx(); ax2.plot(d, [(x - p[0]) / 2 for x in p], "s--", color="#d0182a", label="ΔP por pistón [kW]")
    ax2.set_ylabel("ΔP pérdidas [kW]", color="#d0182a")
    ax[1].axvline(630 - 75 - 345, color="k", ls=":", lw=1); ax[1].text(630 - 75 - 345, 5, " nominal", fontsize=8)
    ax[1].set_xlabel("inserción del pistón Ø40 desde la tapa [mm]"); ax[1].set_ylabel("Δf [kHz]", color="#3b6fd8")
    ax[1].set_title("Sintonizador (carcasa del dee 1)"); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(ROOT / "validation/rf_acoplamiento.png", dpi=140)
print("ok")
