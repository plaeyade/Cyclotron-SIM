"""Techo del búnker del proyecto: dosis en el cuarto de aire acondicionado de encima (v0.28, estudio PRELIMINAR).

Entradas: out/mc_techo/h<altura>/res/t<cm>_s<semilla>.txt (mc_benchmark/bunker_techo, Geant4).
  - corrida de techo grueso: H*(10)-fluencia por capa de 10 cm dentro del techo → curva dosis–espesor;
  - corridas de espesor real: H*(10) en el aire del cuarto de A/A → calibran la curva (la capa incluye el retroceso
    del hormigón de encima).
Criterio (México, Reglamento General de Seguridad Radiológica): público 1/10 del POE (Art. 37: 5 mSv/año) y 1 mSv/año
de promedio para exposición prolongada (Art. 38) → se diseña a 1 mSv/año; zona de radiación > 0.05 mSv/h (Art. 113).
Ocupación del cuarto de A/A: T = 1/40 (ocupación mínima, NCRP 151 tabla 4.1) y 1/20 (almacén) como variante.
Carga: escenarios en µAh/año (todos los blancos). Referencia real: VUB-UZ Brussel ~104 000 µAh/año (VUB24).
Factor SRG: tras 2 m de hormigón el modelo da 0.06–0.08 del SRG (benchmark v0.27) → la variante «escalada al SRG»
multiplica la dosis del modelo por 1/0.07 ≈ 14 (conservadora).
Salida: validation/techo_bunker.{json,png}
"""
import glob, json, math, pathlib, re, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
P_POR_UAH = 1e-6 * 3600 / 1.602176634e-19            # protones por µAh
LIM_ANUAL_uSv = 1000.0                                 # Art. 38 RGSR
ZONA_RAD_uSv_h = 50.0                                  # Art. 113 RGSR
T_OCUP = {"1/40": 1 / 40, "1/20": 1 / 20}
CARGAS = [50_000, 100_000, 150_000]
I_REF_uA = 200.0                                       # 2 × 100 µA (SRG)
B = json.load(open(ROOT / "validation/benchmark_srg.json"))
F_SRG = 1 / np.mean([B["tabla9"][t]["cociente_modelo_SRG_conYsrg"] for t in ("200", "230")])


def leer(f):
    hdr = open(f).readline(); N = int(re.search(r"N (\d+)", hdr).group(1)); tc = float(re.search(r"tc_cm (\S+)", hdr).group(1))
    det, cap = {}, {}
    for l in open(f):
        if l.startswith("# capa "):
            _, _, reg, k, hn, s2n, hg, s2g = l.split(); cap[(reg, int(k))] = (float(hn), float(s2n), float(hg), float(s2g))
        elif not l.startswith("#"):
            x = l.split(); det[x[0]] = (float(x[1]), float(x[2]), float(x[3]), float(x[4]))
    return N, tc, det, cap


def media(lst):
    """lst: [(N, m_n, s2_n, m_g, s2_g)] → (tot, err_rel)"""
    Nt = sum(x[0] for x in lst)
    out = []
    for i in (1, 3):
        m = sum(x[i] * x[0] for x in lst) / Nt; s2 = sum(x[i + 1] for x in lst)
        out.append((m, math.sqrt(max(s2 / Nt - m * m, 0) / Nt)))
    tot = out[0][0] + out[1][0]
    return tot, math.hypot(out[0][1], out[1][1]) / tot if tot > 0 else None, out[1][0] / tot if tot > 0 else None


def uSv_h(pSv_por_n, y):                               # 2 × 100 µA
    return pSv_por_n * y * I_REF_uA * 1e-6 / 1.602176634e-19 * 3600 * 1e-6


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    res = dict(criterio=dict(limite_anual_uSv=LIM_ANUAL_uSv, zona_radiacion_uSv_h=ZONA_RAD_uSv_h, ocupacion=T_OCUP,
                             cargas_uAh_anio=CARGAS, factor_SRG=F_SRG), alturas={})
    for d in sorted(glob.glob(str(ROOT / "out/mc_techo/h*"))):
        geo = json.load(open(pathlib.Path(d) / "geometria.json")); Y = geo["yield_BON21"]
        runs = {}
        for f in glob.glob(d + "/res/t*_s*.txt"):
            N, tc, det, cap = leer(f); runs.setdefault(tc, []).append((N, det, cap))
        if not runs:
            continue
        h = int(pathlib.Path(d).name[1:]); R = dict(puertas=geo["puertas"], interior_mm=geo["inner"], espesores={}, curva={})
        for tc, lst in sorted(runs.items()):
            e = {}
            for k in lst[0][1]:
                tot, er, fg = media([(x[0], *x[1][k]) for x in lst])
                e[k] = dict(pSv_por_n=tot, err_rel=er, frac_gamma=fg, uSv_h_200uA=uSv_h(tot, Y))
            R["espesores"][tc] = dict(semillas=len(lst), detectores=e)
        tmax = max(runs); lst = runs[tmax]
        for reg in ("col_T1", "sala"):
            z, v, ev = [], [], []
            for k in range(int(tmax // 10)):
                if (reg, k) not in lst[0][2]:
                    continue
                tot, er, _ = media([(x[0], *x[2][(reg, k)]) for x in lst])
                if tot > 0 and er is not None and er < 0.5:
                    z.append(10 * k + 10); v.append(tot); ev.append(er)   # dosis «a la salida» de la capa k (espesor 10k+10)
            R["curva"][reg] = dict(espesor_cm=z, pSv_por_n=v, err_rel=ev)
        # calibración: detector de salida / capa equivalente, en los espesores reales corridos
        cal = []
        for tc in runs:
            if tc == tmax:
                continue
            for reg, dn in (("col_T1", "techo_sobre_T1"), ("sala", "techo_promedio_sala")):
                c = R["curva"][reg]
                if tc in c["espesor_cm"]:
                    cal.append(dict(t_cm=tc, region=reg, salida_sobre_capa=R["espesores"][tc]["detectores"][dn]["pSv_por_n"] / c["pSv_por_n"][c["espesor_cm"].index(tc)]))
        # el cociente medido (0.2–1.7, con errores del 30–60 %) confirma que la capa sobreestima la salida; para diseño
        # NO se reduce (factor 1, conservador): la curva por capas es una cota superior
        fcal = 1.0
        R["calibracion"] = dict(casos=cal, mediana_medida=float(np.median([c["salida_sobre_capa"] for c in cal])) if cal else None,
                                factor=fcal, nota="se usa 1 (conservador); la salida real queda ×0.2–0.5 de la curva")
        # diseño: espesor para 1 mSv/año
        dis = []
        c = R["curva"]["col_T1"]; t = np.array(c["espesor_cm"]); v = np.array(c["pSv_por_n"]) * fcal
        a, b = np.polyfit(t, np.log10(v), 1)
        R["TVL_techo_cm"] = -1 / a
        for W in CARGAS:
            for tn, T in T_OCUP.items():
                for var, fs in (("modelo", 1.0), ("escalado_SRG", F_SRG)):
                    # dosis anual = (pSv/n · fs) · n/p · p/µAh · W · T
                    target = LIM_ANUAL_uSv / (Y * P_POR_UAH * W * T * fs * 1e-6)    # pSv/n permitido
                    tc_req = (math.log10(target) - b) / a
                    dis.append(dict(carga_uAh=W, ocupacion=tn, variante=var, espesor_cm=round(tc_req, 1)))
        tz = (math.log10(ZONA_RAD_uSv_h / uSv_h(1.0, Y) / F_SRG) - b) / a
        R["diseno"] = dis; R["espesor_zona_radiacion_cm_escalado_SRG"] = round(tz, 1)
        res["alturas"][h] = R
    json.dump(res, open(ROOT / "validation/techo_bunker.json", "w"), indent=1, ensure_ascii=False, default=float)
    # figura
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for h, R in res["alturas"].items():
        Y = json.load(open(ROOT / f"out/mc_techo/h{h}/geometria.json"))["yield_BON21"]; f = R["calibracion"]["factor"]
        for reg, ls, lab in (("col_T1", "-", "sobre el blanco T1"), ("sala", "--", "promedio de la sala")):
            c = R["curva"][reg]
            ax.semilogy(c["espesor_cm"], [uSv_h(x * f, Y) for x in c["pSv_por_n"]], ls, color="#16a085", label=f"modelo, {lab} (altura {h} mm)")
            ax.semilogy(c["espesor_cm"], [uSv_h(x * f * F_SRG, Y) for x in c["pSv_por_n"]], ls, color="#d0182a", lw=1, label=f"× {F_SRG:.0f} (escalado al SRG), {lab}")
        for tc, E in R["espesores"].items():
            if tc != max(R["espesores"]):
                d = E["detectores"]["techo_sobre_T1"]
                ax.errorbar([tc], [d["uSv_h_200uA"]], [d["uSv_h_200uA"] * d["err_rel"]], fmt="o", color="#16a085", capsize=3)
    for W, col in zip(CARGAS, ("#888", "#555", "#222")):
        r = LIM_ANUAL_uSv / (W / I_REF_uA) / T_OCUP["1/40"]          # tasa a 200 µA que da 1 mSv/año
        ax.axhline(r, color=col, ls=":", lw=1); ax.text(101 if W != 150_000 else 141, r * 1.08, f"1 mSv/año, {W // 1000}k µAh/año, T = 1/40", fontsize=7, color=col)
    ax.axhline(ZONA_RAD_uSv_h, color="#8e44ad", ls="-.", lw=1); ax.text(101, ZONA_RAD_uSv_h / 1.6, "zona de radiación (50 µSv/h, Art. 113)", fontsize=7, color="#8e44ad")
    ax.axvspan(170, 230, color="#3b6fd8", alpha=0.08); ax.text(172, 2e-3, "rango SRG\n1.7–2.3 m", fontsize=7, color="#3b6fd8")
    ax.set_xlim(100, 240); ax.set_ylim(1e-3, 1e3)
    ax.set_xlabel("espesor del techo [cm] (hormigón 2.35 g/cm³)"); ax.set_ylabel("H*(10) en el cuarto de A/A [µSv/h] a 2 × 100 µA")
    ax.set_title("Techo del búnker (preliminar): dosis durante la irradiación"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=6.5, loc="upper right")
    fig.tight_layout(); fig.savefig(ROOT / "validation/techo_bunker.png", dpi=140)
    for h, R in res["alturas"].items():
        print(h, "TVL", round(R["TVL_techo_cm"], 1), "cal", R["calibracion"], "zona", R["espesor_zona_radiacion_cm_escalado_SRG"])
        for x in R["diseno"]:
            print("  ", x)
