"""Región central v0.32: diámetro de la punta de la chimenea (con el bisel de 5 mm del puller), centro fijo a 40.3° del
borde de la punta del dee y pared de la rendija de 0.4 mm. Salidas: validation/region_central_v32.{png,json}"""
import json, math, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; D = ROOT / "out/central_region/puller"
casos = [(2.5, "v32_rc2.5"), (3.0, "v31_ch5"), (4.0, "v32_rc4"), (5.0, "v32_rc5"), (6.5, "v32_rc6.5")]
F = []
for rc, t in casos:
    e = json.load(open(D / t / "resultado.json"))["especies"][0]
    n = sum(e["perdidas_haz"].values()); p = e["perdidas_haz"]["vive"] / n
    lo, hi = e["ventana_fase_aceptada"]; win = hi - lo + 20
    F.append(dict(diametro_mm=2 * rc, tag=t, T=100 * p, err=100 * math.sqrt(p * (1 - p) / n), aceptancia_deg=e["ancho_aceptancia_deg"],
                  dc=100 * p * win / 360, dc_err=100 * math.sqrt(p * (1 - p) / n) * win / 360))
fig, ax = plt.subplots(figsize=(6.6, 4.3))
ax.errorbar([f["diametro_mm"] for f in F], [f["dc"] for f in F], [f["dc_err"] for f in F], fmt="o-", color="#16a085", capsize=3)
for f in F:
    ax.annotate(f"{f['aceptancia_deg']:.0f}°", (f["diametro_mm"], f["dc"]), textcoords="offset points", xytext=(4, 6), fontsize=7)
ax.axvline(6, color="#d0182a", lw=1); ax.text(6.1, 10.8, "elegido Ø6", fontsize=8, color="#d0182a")
ax.axvspan(11, 15, color="#888", alpha=0.12); ax.text(11.2, 15.3, "cuerpo Ø13 ± 2\n(mapa Opera)", fontsize=7, color="#555")
ax.set_xlabel("diámetro de la punta de la chimenea [mm]"); ax.set_ylabel("% de la corriente DC de la fuente a 6 vueltas")
ax.set_title("Chimenea (puller con bisel de 5 mm; etiqueta: fases aceptadas)", fontsize=10); ax.grid(alpha=.3); ax.set_ylim(10, 16)
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_v32.png", dpi=140)
json.dump(dict(casos=F, elegido_diametro_mm=6.0, nota="centro de la chimenea fijo (40.3° antes del borde, r 12.5); pared de la rendija 0.4 mm"),
          open(ROOT / "validation/region_central_v32.json", "w"), indent=1, ensure_ascii=False, default=float)
print([(f["diametro_mm"], round(f["dc"], 2)) for f in F])
