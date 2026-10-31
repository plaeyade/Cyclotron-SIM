"""Efecto de las puertas blindadas de HDPE (18 cm, VUB24) en el búnker del SRG — Geant4 (v0.27).

Compara out/mc_srg/res (sin puertas) con out/mc_srg_puertas/res (puertas en T1 y T5) para los espesores corridos en ambos:
  - dosis H*(10) en los detectores interiores y tras el muro (cociente sin/con puertas);
  - traza ponderada de neutrones por capa de 10 cm del muro frontal (∝ fluencia integrada en la capa), como
    indicador de la activación del hormigón; IBA calculó ÷2.3 a 0–10 cm con puertas (VUB24) y el SRG §2.1.8.8 da ÷10
    en el flujo de neutrones.
Salida: validation/benchmark_puertas.{json,png}
"""
import glob, json, math, pathlib, re
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
CASOS = {"sin": ROOT / "out/mc_srg/res", "con": ROOT / "out/mc_srg_puertas/res"}
REF = dict(IBA_activacion_0_10cm=2.3, SRG_flujo_n=10.0)


def leer(f):
    N = int(re.search(r"N (\d+)", open(f).readline()).group(1))
    det, capa = {}, {}
    for l in open(f):
        if l.startswith("# capa "):
            _, _, k, n, w = l.split(); capa[int(k)] = float(w)
        elif not l.startswith("#"):
            x = l.split(); det[x[0]] = (float(x[1]) + float(x[3]), float(x[1]))
    return N, det, capa


def junta(carpeta, t):
    fs = sorted(glob.glob(str(carpeta / f"t{t}_s*.txt")))
    if not fs:
        return None
    L = [leer(f) for f in fs]
    out = dict(semillas=len(L), det={}, capa={})
    for d in L[0][1]:
        v = [x[1][d][0] for x in L]; vn = [x[1][d][1] for x in L]
        out["det"][d] = dict(tot=float(np.mean(v)), n=float(np.mean(vn)),
                             disp=float(np.std(v) / np.mean(v) / math.sqrt(len(v))) if len(v) > 1 and np.mean(v) > 0 else None)
    for k in L[0][2]:
        v = [x[2].get(k, 0.0) for x in L]
        out["capa"][k] = dict(w=float(np.mean(v)), disp=float(np.std(v) / np.mean(v) / math.sqrt(len(v))) if len(v) > 1 and np.mean(v) > 0 else None)
    return out


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    ts = sorted({int(re.search(r"t(\d+)_s", f).group(1)) for f in glob.glob(str(CASOS["con"] / "t*_s*.txt"))})
    res = dict(referencias=REF, espesores={})
    for t in ts:
        a, b = junta(CASOS["sin"], t), junta(CASOS["con"], t)
        if not (a and b):
            continue
        r = dict(semillas=dict(sin=a["semillas"], con=b["semillas"]), dosis={}, capas={})
        for d in a["det"]:
            A, B = a["det"][d], b["det"][d]
            r["dosis"][d] = dict(sin_pSv_por_n=A["tot"], con_pSv_por_n=B["tot"],
                                 reduccion_total=A["tot"] / B["tot"] if B["tot"] > 0 else None,
                                 reduccion_neutrones=A["n"] / B["n"] if B["n"] > 0 else None)
        for k in sorted(a["capa"]):
            if k in b["capa"] and b["capa"][k]["w"] > 0:
                r["capas"][k] = dict(prof_cm=[10 * k, 10 * k + 10], sin=a["capa"][k]["w"], con=b["capa"][k]["w"],
                                     reduccion=a["capa"][k]["w"] / b["capa"][k]["w"])
        res["espesores"][t] = r
    # perfil en profundidad (con puertas, muro más grueso) frente a las varillas de hormigón medidas en VUB24
    tmax = max(res["espesores"]); cp = res["espesores"][tmax]["capas"]
    kk = [k for k in sorted(cp) if k < 9]; z = np.array([10 * k + 5 for k in kk]); w = np.array([cp[k]["con"] for k in kk])
    tramo = lambda a, b: float(np.mean(np.interp(np.linspace(a, b, 31), z, w)))
    perfil = [tramo(a, a + 15) for a in (0, 15, 30, 45)]
    act = json.load(open(ROOT / "validation/activacion_vub_kiube.json"))["medida"]
    med = {f"{u}_{n}": [round(x / [y for y in v if y][0], 3) if x else None for x in v]
           for u, d in act.items() if u.startswith("ubicacion") for n, v in d.items()}
    res["perfil_profundidad"] = dict(tramos_cm=["0-15", "15-30", "30-45", "45-60"], muro_cm=tmax,
                                     modelo_con_puertas_rel=[round(p / perfil[0], 3) for p in perfil], medido_VUB24_rel=med,
                                     nota="modelo: fluencia de neutrones (todas las energías); medida: Eu-152 y Co-60 (captura térmica)")
    json.dump(res, open(ROOT / "validation/benchmark_puertas.json", "w"), indent=1, ensure_ascii=False, default=float)
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for t, r in res["espesores"].items():
        ks = [k for k in sorted(r["capas"]) if k < min(int(t) // 10, 10)]; z = [10 * k + 5 for k in ks]   # >1 m: estadística pobre
        ax[0].semilogy(z, [r["capas"][k]["reduccion"] for k in ks], "o-", label=f"muro {t} cm")
    ax[0].axhline(REF["IBA_activacion_0_10cm"], color="#d0182a", ls="--", lw=1, label="IBA: activación 0–10 cm ÷2.3 (VUB24)")
    ax[0].axhline(REF["SRG_flujo_n"], color="#888", ls=":", lw=1, label="SRG §2.1.8.8: flujo n ÷10")
    ax[0].set_xlabel("profundidad en el muro frontal [cm]"); ax[0].set_ylabel("sin puertas / con puertas")
    ax[0].set_title("Fluencia de neutrones en el muro frontal (primer metro)"); ax[0].set_ylim(1, 20); ax[0].grid(alpha=.3, which="both"); ax[0].legend(fontsize=7.5)
    dets = ["frente_0deg_interior", "a_1m_del_blanco_0deg", "lateral_90deg_interior", "opuesto_180deg_interior", "exterior_tras_muro"]
    w = 0.8 / max(len(res["espesores"]), 1)
    for i, (t, r) in enumerate(res["espesores"].items()):
        ax[1].bar(np.arange(len(dets)) + i * w, [r["dosis"][d]["reduccion_total"] or 0 for d in dets], w, label=f"muro {t} cm")
    ax[1].set_xticks(np.arange(len(dets)) + 0.4 - w / 2)
    ax[1].set_xticklabels(["0° junto\nal muro", "1 m del\nblanco", "90°", "180°", "tras el\nmuro"], fontsize=8)
    ax[1].set_yscale("log"); ax[1].axhline(1, color="k", lw=0.6); ax[1].set_ylabel("H*(10) sin / con puertas")
    ax[1].set_title("Reducción de dosis con las puertas de HDPE de 18 cm"); ax[1].grid(alpha=.3, axis="y", which="both"); ax[1].legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(ROOT / "validation/benchmark_puertas.png", dpi=140)
    for t, r in res["espesores"].items():
        print(t, {d: round(v["reduccion_total"], 2) if v["reduccion_total"] else None for d, v in r["dosis"].items()})
        print("   capas", {k: round(v["reduccion"], 2) for k, v in r["capas"].items()})
