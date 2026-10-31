"""Fuente → blanco (v0.35): balance de partículas, nivel de campo, pérdidas verticales y mancha en el blanco.
Lee out/fuente_blanco/ (fem/fuente_a_blanco.py). Salidas: validation/fuente_blanco.{png,json}"""
import json, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; D = ROOT / "out/fuente_blanco"
GAS = 0.15                                                  # pérdida por gas residual (vacuum: calibrada con TUD03, 81 %)
VENTANA = 120.0                                             # fases de emisión simuladas (−80…+40°) de 360°

rc = np.load(D / "rc_k1.0103.npz", allow_pickle=True)
A = json.load(open(D / "acel_p0_final.json")); X = json.load(open(D / "ext_p0_final.json"))
n0 = int(rc["n"]); m_rc = rc["motivo"]
n_rc = int(np.isin(m_rc, ["vive", "r"]).sum())
mot = np.array(A["motivo"]); n_lam = int((mot == "lamina").sum())
f_col, f_win = X["frac_colimador"], X["frac_ventana"]
# fracciones de la corriente DC de la fuente
dc = VENTANA / 360
etapas = [("emitido en la ventana\n−80…+40° (DC × 1/3)", dc),
          ("sale de la región\ncentral (r 126 mm)", dc * n_rc / n0),
          ("llega a la lámina\n(sin gas)", dc * n_lam / n0),
          ("llega a la lámina\n(con gas, −15 %)", dc * n_lam / n0 * (1 - GAS)),
          ("dentro de la ventana\ndel blanco Ø20", dc * n_lam / n0 * (1 - GAS) * f_win),
          ("pasa el colimador\nØ9", dc * n_lam / n0 * (1 - GAS) * f_col)]

fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
a = ax[0, 0]
y = np.arange(len(etapas))[::-1]
a.barh(y, [100 * v for _, v in etapas], color=["#bbb", "#8e44ad", "#3b6fd8", "#16a085", "#e67e22", "#d0182a"])
for yy, (lab, v) in zip(y, etapas):
    a.text(100 * v + 0.4, yy, f"{100 * v:.1f} %", va="center", fontsize=8)
a.set_yticks(y); a.set_yticklabels([l for l, _ in etapas], fontsize=7.5)
a.set_xlabel("% de la corriente DC extraída de la fuente"); a.set_xlim(0, 40)
a.set_title(f"Balance fuente → blanco (puerto 1, nivel ×1.0103, {n0} partículas)", fontsize=9.5)

a = ax[0, 1]
ks, fr = [], []
for k in ("1.0099", "1.0101", "1.0103", "1.0105", "1.0107"):
    f = D / f"acel_p0_n800_k{k}.json"
    if f.exists():
        R = json.load(open(f)); ks.append(float(k)); fr.append(100 * R["motivos"].get("lamina", 0) / R["n"])
a.plot([100 * (k - 1) for k in ks], fr, "o-", color="#3b6fd8")
a.set_xlabel("nivel de campo sobre el mapa de 115 kA·v [%]"); a.set_ylabel("% de la salida de la región central que llega a la lámina")
a.set_title("Ajuste de la corriente del imán (como en operación)", fontsize=9.5); a.grid(alpha=.3)

a = ax[1, 0]
Ez = [h[-1][1] for h, m in zip(A["hist"], A["motivo"]) if m == "z" and h]
a.hist(Ez, bins=np.arange(0, 18.6, 1.0), color="#d0182a", alpha=.8, label=f"pérdidas verticales ({len(Ez)})")
a.set_xlabel("energía de la pérdida [MeV]"); a.set_ylabel("partículas")
b = a.twinx()
E = json.load(open(ROOT / "out/orbitas/equilibrio_v34.json"))
T = [o["T_MeV"] for o in E]; d = [o["nu_r"] - 2 * o["nu_z"] for o in E]
b.plot(T, d, "k.-", label="FE v0.34: νr − 2νz")
K = json.load(open(ROOT / "validation/isocronia.json"))[0]
b.plot(K["T_MeV"], np.array(K["nu_r"]) - 2 * np.array(K["nu_z"]), "--", color="#3b6fd8", label="KEY22 medido: νr − 2νz")
b.axhline(0, color="#888", lw=.8); b.set_ylim(-1.1, 0.4); b.set_ylabel("νr − 2νz (resonancia de Walkinshaw en 0)")
h1, l1 = a.get_legend_handles_labels(); h2, l2 = b.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="center left")
a.set_title("Pérdidas verticales: cruce de νr = 2νz a 14–17 MeV", fontsize=9.5)

a = ax[1, 1]
p = X["puntos"]; ok = np.isfinite(p["x"])
a.scatter(p["x"], p["z"], s=8, c=p["T"], cmap="viridis")
for rr, c, lab in ((4.5, "#d0182a", "colimador Ø9"), (10, "#e67e22", "ventana Ø20")):
    t = np.linspace(0, 2 * np.pi, 200); a.plot(rr * np.cos(t), rr * np.sin(t), color=c, label=lab)
a.set_aspect("equal"); a.set_xlim(-35, 35); a.set_ylim(-35, 35); a.grid(alpha=.3); a.legend(fontsize=7.5)
a.set_xlabel("horizontal en el blanco [mm]"); a.set_ylabel("vertical [mm]")
a.set_title(f"Mancha en el blanco: σx {X['sigma_x_mm']:.1f} mm, σz {X['sigma_z_mm']:.1f} mm; T {X['T_MeV'][0]:.2f} ± {X['T_MeV'][1]:.2f} MeV",
            fontsize=9.5)
fig.tight_layout(); fig.savefig(ROOT / "validation/fuente_blanco.png", dpi=130)
res = dict(n_emitidas=n0, n_salida_rc=n_rc, n_lamina=n_lam, motivos_acel=A["motivos"], perdidas_rc={x: int((m_rc == x).sum()) for x in set(m_rc)},
           fraccion_DC=[dict(etapa=l.replace("\n", " "), pct=100 * v) for l, v in etapas], nivel_scan=dict(nivel=ks, pct_lamina=fr),
           extraccion={k: v for k, v in X.items() if k not in ("puntos", "barrido")},
           perdidas_z_14_18_MeV=int(sum(1 for e in Ez if e >= 14)), perdidas_z=len(Ez))
json.dump(res, open(ROOT / "validation/fuente_blanco.json", "w"), indent=1, default=float, ensure_ascii=False)
print(json.dumps(res["fraccion_DC"], ensure_ascii=False))
