"""Menisco del plasma (v0.30): forma autoconsistente por nivel de voltaje y transmisión de la región central.
Entradas: out/central_region/menisco/menisco_jp*.json y transmision_jp*.json
Salidas: validation/menisco.{png,json}"""
import glob, json, math, pathlib, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "out/central_region/menisco"
W, H, RC, PD = 0.6, 6.0, 3.0, 2.6

def sm(c, z, u=0.0):
    z2 = np.asarray(z) ** 2
    return c[0] + c[1] * z2 + c[2] * z2 ** 2 + c[3] * z2 ** 3 + (c[4] if len(c) > 4 else 0) * (2 * np.asarray(u) / W) ** 2

M = json.load(open(D / "menisco_jp1440.json"))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
z = np.linspace(-H / 2, H / 2, 121); u = np.linspace(-W / 2, W / 2, 61)
cols = ["#16a085", "#2e86c1", "#8e44ad", "#d0182a"]
for r, col in zip(sorted(M["niveles"], key=lambda r: -r["cos"]), cols):
    lab = f"V = {r['V'] / 1e3:.1f} kV (cos φ = {r['cos']})"
    ax[0].plot(sm(r["coef"], z), z, color=col, label=lab)
    ax[1].plot(u, sm(r["coef"], 0.0, u), color=col, label=lab)
for a in ax[:2]:
    a.axvline(PD, color="#888", ls="--", lw=1); a.axvline(RC, color="k", lw=1.2)
ax[0].text(PD + 0.01, -2.9, "cara plana\n(v0.29)", fontsize=7, color="#666"); ax[0].text(RC + 0.01, 2.6, "boca de la\nrendija", fontsize=7)
ax[0].set_xlabel("s: profundidad desde el eje de la chimenea [mm]"); ax[0].set_ylabel("z [mm]"); ax[0].set_xlim(2.2, 3.05)
ax[0].set_title("Menisco en el plano vertical (u = 0)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
ax[1].set_xlabel("s [mm]"); ax[1].set_ylabel("u: a lo ancho de la rendija [mm]")
ax[1].cla()
for r, col in zip(sorted(M["niveles"], key=lambda r: -r["cos"]), cols):
    ax[1].plot(sm(r["coef"], 0.0, u), u, color=col)
ax[1].axvline(PD, color="#888", ls="--", lw=1); ax[1].axvline(RC, color="k", lw=1.2)
ax[1].axhline(W / 2, color="k", lw=2); ax[1].axhline(-W / 2, color="k", lw=2)
ax[1].set_xlabel("s [mm]"); ax[1].set_ylabel("u [mm]"); ax[1].set_xlim(2.2, 3.05); ax[1].set_ylim(-0.4, 0.4)
ax[1].set_title("Menisco a lo ancho (z = 0): cóncavo, paredes en ±0.3"); ax[1].grid(alpha=.3)
# transmisión
casos = []
for f in sorted(glob.glob(str(D / "transmision_jp*.json"))):
    R = json.load(open(f)); jp = R["jp_A_m2"]
    casos.append(dict(etiqueta=("plasma plano (control)" if R["plano"] else f"menisco, j_p {jp / 10:.0f} mA/cm²"),
                      jp=jp, plano=R["plano"], T=R["transmision_pct"], err=R["error_pct"], perdidas=R["perdidas_pct"]))
casos.sort(key=lambda c: (not c["plano"], c["jp"]))
x = np.arange(len(casos))
ax[2].bar(x, [c["T"] for c in casos], yerr=[c["err"] for c in casos], color=["#888" if c["plano"] else "#16a085" for c in casos], capsize=4)
for i, c in enumerate(casos):
    ax[2].text(i, 2, f"z {c['perdidas']['z']:.0f} %\nmetal {c['perdidas']['metal']:.0f} %", ha="center", fontsize=7, color="w")
ax[2].set_xticks(x); ax[2].set_xticklabels([c["etiqueta"] for c in casos], fontsize=7, rotation=12)
ax[2].set_ylabel("transmisión 6 vueltas [%]"); ax[2].set_title("Región central: −80…+40° RF, 400 partículas"); ax[2].grid(alpha=.3, axis="y")
fig.tight_layout(); fig.savefig(ROOT / "validation/menisco.png", dpi=130)
json.dump(dict(niveles=[dict(cos=r["cos"], V=r["V"], coef=r["coef"], s_centro=r["s_nodos"][0], a_u=r["coef"][4],
                             iteraciones=len(r["historia"]), salida=r["salida"]) for r in M["niveles"]], transmision=casos,
               jp_ref="Kuo & Hendry, Cyclotrons 1992: 144 mA/cm² (2.8 mA DC equivalentes, rendija 0.6 × 4 mm)"),
          open(ROOT / "validation/menisco.json", "w"), indent=1, ensure_ascii=False, default=float)
print([(c["etiqueta"], round(c["T"], 1)) for c in casos])
