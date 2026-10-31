"""Región central v0.29: giro de la rendija de la chimenea, forma de la nariz del dee (radio de punta, ángulo) y semiabertura
de las puntas, con el cuello entre puntas (TUD03T) — preset kiube, malla fina, 6 vueltas, 400 partículas.
Entradas: out/central_region/puller/v29*/resultado.json y out/central_region/diag_z/diag.json
Salidas: validation/region_central_v29.{png,json}"""
import json, math, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "out/central_region/puller"


def fila(tag):
    e = json.load(open(D / tag / "resultado.json"))["especies"][0]
    n = sum(e["perdidas_haz"].values()); p = e["perdidas_haz"]["vive"] / n
    return dict(tag=tag, T=100 * p, err=100 * math.sqrt(p * (1 - p) / n), z=100 * e["perdidas_haz"]["z"] / n,
                metal=100 * e["perdidas_haz"]["metal"] / n, acept=e["ancho_aceptancia_deg"])


base = fila("v29_base")
series = {
    "giro de la rendija [°] (+ = hacia fuera)": [(-20, "v29_rot-20"), (-10, "v29_rot-10"), (-5, "v29_rot-5"), (0, "v29_base"),
                                                 (5, "v29_rot+5"), (10, "v29_rot+10"), (20, "v29_rot+20")],
    "radio de la punta del dee [mm]": [(6, "v29_nose6"), (8, "v29_base"), (10, "v29_nose10")],
    "ángulo de las puntas [°]": [(44, "v29_ang44"), (50, "v29_base"), (56, "v29_ang56")],
    "semiabertura de las puntas [mm]": [(7, "v29_ap7"), (8, "v29_base"), (9, "v29_ap9"), (10, "v29_ap10")],
}
out = dict(base=base, series={}, sin_cuello=fila("v29_sinpuente"), modelo_corregido=fila("v29b_base"))
fig, ax = plt.subplots(1, 5, figsize=(19, 3.9))
for a, (lab, pts) in zip(ax, series.items()):
    F = [dict(x=x, **fila(t)) for x, t in pts]; out["series"][lab] = F
    a.errorbar([f["x"] for f in F], [f["T"] for f in F], [f["err"] for f in F], fmt="o-", color="#16a085", capsize=3)
    a.axhspan(base["T"] - base["err"], base["T"] + base["err"], color="#16a085", alpha=0.12)
    a.set_xlabel(lab); a.set_ylim(0, 45); a.grid(alpha=.3)
ax[0].set_ylabel("transmisión 6 vueltas [%]")
dz = json.load(open(ROOT / "out/central_region/diag_z/diag.json"))
v = [f["vuelta"] for f in dz["filas"] if f["motivo"] == "z"]
ax[4].hist(v, bins=24, color="#d0182a"); ax[4].set_xlabel("vuelta en que se pierde por |z|"); ax[4].set_title(
    f"{100 * dz['perdidas_z'] / dz['n']:.0f} % del haz, mediana {np.median(v):.2f} vueltas", fontsize=9)
fig.suptitle("Región central v0.29 (con cuello entre puntas): ninguna variable geométrica sale de la meseta; "
             "la pérdida es vertical en el primer cuarto de vuelta (óptica de extracción de la rendija)", fontsize=10)
fig.tight_layout(); fig.savefig(ROOT / "validation/region_central_v29.png", dpi=130)
out["diag_z"] = {k: v for k, v in dz.items() if k != "filas"}
json.dump(out, open(ROOT / "validation/region_central_v29.json", "w"), indent=1, ensure_ascii=False, default=float)
print("ok", base, out["modelo_corregido"])
