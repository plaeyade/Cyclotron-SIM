"""Calibración del circuito magnético (v0.18): resume los casos FE de out/fem_mag/calib y variantes frente a KEY22.

Métricas en la banda r = 0.38–0.47 m (como KEY22 Fig. 3): ⟨B⟩ medio y flutter medio; además B_colina/B_valle en 0.45 m,
el reparto de Ampère por tramos (si el caso se corrió con --diag) y el chequeo BATAN13: B_valle / (μ0·NI/g_valle).
Uso: python validation/calibrar_iman.py  → validation/calibracion_iman.json + .png
"""
import glob, json, math, pathlib, re
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
MU0 = 4e-7 * math.pi
d = np.genfromtxt(ROOT / "validation/kiube_Bavg_flutter_KEY22.csv", delimiter=",", skip_header=4)
banda = (d[:, 0] >= 0.38) & (d[:, 0] <= 0.47)
KEY = dict(B_banda_T=float(np.nanmean(d[banda, 1]) / 1e4), F_banda=float(np.nanmean(d[banda, 2])),
           B_centro_rel=float(np.nanmax(d[:12, 1]) / np.interp(0.45, d[~np.isnan(d[:, 1]), 0], d[~np.isnan(d[:, 1]), 1])))
S09 = dict(B_colina_T=2.15, B_valle_T=0.45)
BATAN = 0.769 / (MU0 * 88024 / 0.120)                    # 0.84


def caso(f):
    R = json.load(open(f)); r = np.array(R["r_m"]); B = np.array(R["Bavg_T"]); F = np.array(R["flutter"])
    m = (r >= 0.38) & (r <= 0.47)
    out = dict(NI_total=R["NI_total"], B_banda_T=float(np.nanmean(B[m])), F_banda=float(np.nanmean(F[m])),
               B_colina_045=float(np.interp(0.45, r, R["B_colina_T"])), B_valle_045=float(np.interp(0.45, r, R["B_valle_T"])),
               B_centro_rel=float(np.interp(0.02, r, B) / np.interp(0.45, r, B)))
    return out


def recolectar():
    filas = []
    for res in sorted(glob.glob(str(ROOT / "out/fem_mag/calib/*/resumen.json"))):
        carpeta = pathlib.Path(res).parent
        meta = {(round(x["NI_total"])): x for x in json.load(open(res))}
        for f in sorted(carpeta.glob("resultado_NI*.json")):
            c = caso(f)
            m = meta.get(round(c["NI_total"]), {})
            c.update(caso_id=carpeta.name, Js=m.get("Js"), ampere=m.get("ampere"))
            filas.append(c)
    return filas


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    filas = recolectar()
    gv = 0.220                                              # gap de valle completo [m] (2 × valley_floor_z) — por caso abajo
    for c in filas:
        vz = re.search(r"valle(\d+)", c["caso_id"])
        vz2 = re.search(r"_v(\d+)", c["caso_id"])
        g = 2 * (float(vz.group(1)) if vz else float(vz2.group(1)) if vz2 else 100.0 if c["caso_id"].startswith("cal20")
                 else 95.0 if c["caso_id"].startswith("cal") else 110.0) / 1e3
        c["eficiencia_valle_BATAN"] = c["B_valle_045"] / (MU0 * c["NI_total"] / g)
        if c.get("ampere"):
            tot = sum(t["Hdl_hierro"] + t["Hdl_aire"] for t in c["ampere"])
            c["fraccion_hierro"] = sum(t["Hdl_hierro"] for t in c["ampere"]) / tot
    # sensibilidad local d ln⟨B⟩ / d ln NI en el punto de trabajo, frente al C11 de IBA (TOSCA, CY11):
    # +160 G por 10 A a 186 A con B0 ≈ 1.408 T (42.877 MHz, h = 2) → (0.016/1.408)/(10/186) = 0.21
    cy11 = (0.0160 / 1.408) / (10 / 186)
    par = [c for c in filas if c["caso_id"] == "cal_p181_v95"]
    exp_fe = None
    if len(par) >= 2:
        a, b = sorted(par, key=lambda c: c["NI_total"])[:2]
        exp_fe = math.log(b["B_banda_T"] / a["B_banda_T"]) / math.log(b["NI_total"] / a["NI_total"])
        print(f"d ln⟨B⟩/d ln NI: FE calibrado {exp_fe:.3f}  ·  C11 de IBA (TOSCA) {cy11:.3f}")
    out = dict(KEY22=KEY, S09=S09, BATAN13_eficiencia_valle=BATAN,
               sensibilidad_dlnB_dlnNI=dict(FE_calibrado=exp_fe, CY11_TOSCA=cy11), casos=filas)
    json.dump(out, open(ROOT / "validation/calibracion_iman.json", "w"), indent=1, ensure_ascii=False)
    for c in filas:
        print(f"{c['caso_id']:>22} NI {c['NI_total'] / 1e3:5.0f}k Js {c['Js']}  ⟨B⟩ {c['B_banda_T']:.4f}  F {c['F_banda']:.3f}  "
              f"Bh {c['B_colina_045']:.2f} Bv {c['B_valle_045']:.3f}  ef.valle {c['eficiencia_valle_BATAN']:.2f}  "
              f"hierro {c.get('fraccion_hierro', float('nan')):.2f}")
    print("KEY22", KEY, "BATAN", round(BATAN, 3))
    # figura: KEY22 frente al FE antes (1010 genérica, 170 kA·v) y después (calibrado, 116 kA·v)
    sel = [("v0.7–v0.17: 1010 genérica, tapa 166, valle 110, 170 kA·v", "base/resultado_NI170k.json", "#888888"),
           ("v0.7–v0.17 a 105 kA·v (la potencia publicada)", "base/resultado_NI105k.json", "#d0182a"),
           ("v0.18: acero J_s 2.12 T, tapa 181, valle 95, 116 kA·v", "cal_final/resultado_NI116k.json", "#16a085"),
           ("v0.21: + valle 100 (PAT16) y cuña del corrector, 115 kA·v = 19 kW", "cal20_final/resultado_NI115k.json", "#8e44ad")]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax[0].plot(d[:, 0], d[:, 1] / 1e4, "o", ms=2, color="#3b6fd8", label="IBA KEY22 (KIUBE)")
    ax[1].plot(d[:, 0], d[:, 2], "o", ms=2, color="#3b6fd8", label="IBA KEY22 (KIUBE)")
    for lab, f, col in sel:
        f = ROOT / "out/fem_mag/calib" / f
        if f.exists():
            R = json.load(open(f))
            ax[0].plot(R["r_m"], R["Bavg_T"], color=col, lw=1.5, label=lab)
            ax[1].plot(R["r_m"], R["flutter"], color=col, lw=1.5, label=lab)
    ax[0].set_ylabel("⟨B⟩(r) [T]"); ax[1].set_ylabel("flutter F")
    for x in ax:
        x.set_xlabel("r [m]"); x.grid(alpha=.3); x.legend(fontsize=7.5); x.set_xlim(0, 0.56)
    ax[0].set_ylim(1.05, 1.45); ax[1].set_ylim(0, 0.5)
    ax[0].set_title("Campo medio"); ax[1].set_title("Flutter")
    fig.tight_layout(); fig.savefig(ROOT / "validation/calibracion_iman.png", dpi=140)
