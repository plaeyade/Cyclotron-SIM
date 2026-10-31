"""Figura de la comparación Molflow 2.11.1 frente al TPMC (lee out/molflow_vs_tpmc/)."""
import json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube.tpmc import TPMC, tube_mesh, sticking_for_speed            # noqa: E402
W = ROOT / "out/molflow_vs_tpmc"; C = json.load(open(W / "comparacion.json")); Z = np.load(W / "tutorial_curvas.npz")
# curva TPMC del tutorial (misma receta que validation/molflow_tutorial_pumpdown.py)
m, c1, c2 = tube_mesh(R=1.0, L=1.0, n=96)
s = np.zeros(len(m.faces))
for cap in (c1, c2):
    s[cap] = sticking_for_speed(1.0, m.area_faces[cap].sum(), 28.0)
dt = 0.05e-3; edges = np.arange(0, 30e-3 + dt / 2, dt)
res = TPMC(m, s, seed=5).run_tof(c1.astype(float) * m.area_faces, ~(c1 | c2), edges, 28.0, n_particles=200000, batch=50000)
h = res["H"] / 2e5; t = 0.5 * (edges[1:] + edges[:-1]); kT = 1.380649e-23 * 293.15
Nd = np.where((t >= 0.01) & (t < 0.02), 1e-4 * 0.1 / kT, 0.0)
Pt = np.convolve(Nd, h)[:len(t)] / (m.area_faces[~(c1 | c2)].sum() * 1e-4) * math.sqrt(2 * math.pi * 28 * 1.66054e-27 * kT) / 100
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4), gridspec_kw=dict(width_ratios=[1, 1.4, 1]))
A = C["A_clausing"]; x = np.arange(len(A))
for k, (key, col, lab) in enumerate((("Clausing", "#444", "Clausing–Cole"), ("Molflow", "#c0392b", "Molflow 2.11.1"),
                                     ("TPMC", "#1f6fb2", "TPMC propio"))):
    ax[0].bar(x + (k - 1) * 0.26, [r[key] for r in A], 0.26, color=col, label=lab)
ax[0].set_xticks(x, [f"L/R = {r['L_R']:g}" for r in A]); ax[0].set_ylabel("transmisión"); ax[0].legend(fontsize=8)
ax[0].set_title("A. Tubos (Clausing)")
ax[1].plot(Z["t"] * 1e3, Z["P_molflow"], color="#c0392b", lw=1.6, label="Molflow 2.11.1 (momentos)")
ax[1].plot(t * 1e3, Pt, color="#1f6fb2", lw=1.0, ls="--", label="TPMC propio (tiempo de vuelo)")
ax[1].plot(Z["t"] * 1e3, Z["P_analitico"], color="k", lw=0.8, ls=":", label="analítico V/S")
ax[1].set_xlabel("t [ms]"); ax[1].set_ylabel("P [mbar]"); ax[1].legend(fontsize=8); ax[1].set_xlim(0, 30)
ax[1].set_title("B. Pump-down del tutorial de Molflow (inyección 10–20 ms)")
K = C["C_kiube"]; lab = ["H₂O de metales, 10 h", "H₂ de la fuente (2 sccm)"]
v_m = [K["H2O_10h"]["P_haz_Molflow"], K["H2_fuente"]["P_haz_Molflow"]]; v_t = [K["H2O_10h"]["P_haz_TPMC"], K["H2_fuente"]["P_haz_TPMC"]]
xx = np.arange(2)
ax[2].bar(xx - 0.18, v_m, 0.36, color="#c0392b", label="Molflow 2.11.1"); ax[2].bar(xx + 0.18, v_t, 0.36, color="#1f6fb2", label="TPMC")
for i in range(2):
    ax[2].text(xx[i], max(v_m[i], v_t[i]) * 1.4, f"{100 * (v_m[i] / v_t[i] - 1):+.1f} %", ha="center", fontsize=9)
ax[2].set_yscale("log"); ax[2].set_xticks(xx, lab, fontsize=8); ax[2].set_ylabel("P en la zona del haz [mbar]")
ax[2].set_title(f"C. Cámara del KIUBE ({K['facetas']:,} facetas)"); ax[2].legend(fontsize=8); ax[2].set_ylim(1e-9, 1e-4)
for a_ in ax:
    a_.grid(alpha=.3)
fig.tight_layout(); fig.savefig(ROOT / "validation/molflow_vs_tpmc.png", dpi=130)
print("ok")
