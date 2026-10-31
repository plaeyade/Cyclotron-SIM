"""Figura del equivalente v0.35 de F70: H*(10) en el plano medio (KIUBE 180, 2 × 90 µA, puertas HDPE) y perfiles a
través de los muros. Lee out/mc_f70v035/<caso>/mapa_h10.npz y resultado_mapa.json → validation/f70_v035_mapa.png"""
import json, pathlib, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ROOT = pathlib.Path(__file__).resolve().parents[1]
d = ROOT / (sys.argv[1] if len(sys.argv) > 1 else "out/mc_f70v035/h2500_t1450")
Z = np.load(d / "mapa_h10.npz"); R = json.load(open(d / "resultado_mapa.json"))
k = float(Z["k_uSv_h_por_pSv_n"]); H = Z["H_media_pSv_por_n"] * k; S = Z["H_std_pSv_por_n"] * k
x, y, z = Z["x_mm"] / 1000, Z["y_mm"] / 1000, Z["z_mm"] / 1000
iz = int(np.argmin(np.abs(z)))
b = R["bloque_hormigon"]; hx, hy, w = b["interior_mm"][0] / 2000, b["interior_mm"][1] / 2000, b["muro_mm"] / 1000
fig = plt.figure(figsize=(13, 5.6))
a = fig.add_subplot(1, 2, 1)
M = np.ma.masked_less_equal(H[:, :, iz].T, 0)
im = a.pcolormesh(x, y, M, norm=LogNorm(vmin=max(M.min(), 1e-3), vmax=M.max()), cmap="inferno", shading="nearest")
for (u, v) in ((hx, hy), (hx + w, hy + w)):
    a.plot([-u, u, u, -u, -u], [-v, -v, v, v, -v], color="#9ad", lw=0.8)
for s in R["fuente"]["blancos"]:
    pass
a.set_aspect("equal"); a.set_xlabel("x [m] (oeste → este)"); a.set_ylabel("y [m] (sur → norte)")
a.set_title(f"H*(10) en el plano medio [µSv/h] — KIUBE 180 (2×90 µA), puertas HDPE\n{R['historias_totales']:.0f} historias, {R['corridas']} corridas", fontsize=9)
fig.colorbar(im, ax=a, fraction=0.04)
a = fig.add_subplot(1, 2, 2)
cols = {"oeste": "#3b6fd8", "este": "#16a085", "sur": "#e67e22", "norte": "#d0182a"}
for nombre in ("oeste", "este", "sur", "norte"):
    P = R["caras_plano_medio"][nombre]["perfil"]
    dd = np.array([f["prof_cm"] for f in P]) / 100; v = np.array([f["cara"] for f in P]); e = np.array([f["cara_se"] for f in P])
    ok = (v > 0) & (e < v)
    tvl = R["caras_plano_medio"][nombre]["tvl_cara_cm"]
    a.errorbar(dd[ok], v[ok], e[ok], fmt="o-", ms=3, color=cols[nombre], capsize=2,
               label=f"{nombre} (TVL {tvl:.0f} cm)" if tvl else nombre)
a.axvspan(0, w, color="#bbb", alpha=0.3); a.text(w / 2, 3e-4, "hormigón 2.11 m", ha="center", fontsize=8)
a.set_yscale("log"); a.set_xlabel("profundidad desde la cara interior del muro [m]"); a.set_ylabel("H*(10) [µSv/h]")
a.axhline(0.5, color="#555", ls=":", lw=0.8); a.text(-0.4, 0.6, "0.5 µSv/h", fontsize=7)
a.set_title("H*(10) medio sobre la cara de cada muro, por capa de 10 cm (± error entre corridas)", fontsize=9); a.grid(alpha=.3, which="both"); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(ROOT / (sys.argv[2] if len(sys.argv) > 2 else "validation/f70_v035_mapa.png"), dpi=130)
print("ok")
