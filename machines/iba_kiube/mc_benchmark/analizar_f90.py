"""Análisis F90: validación contra IBA SRG Rev G tabla 9 y doble haz (2×150 µA) sin puertas blindadas de blanco.
Uso: python analizar_f90.py <tag> [rho_corrida] — lee out/mc_f90/res/<tag>_<cara>_s*_malla.bin
Normalización: 300 µA totales, rendimiento IBA 6.9e-3 n/p (SRG §2.1.8.4) y, como referencia, BON21 4.93e-3 n/p.
Salida: out/mc_f90/<tag>.json con perfiles, puntos IBA (200/220/230 cm de concreto 2.35 por masa) y espesores
requeridos para tabique 2.5 (escalado por masa desde la densidad de la corrida)."""
import glob, json, math, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import leer_iba
import analizar_iba as A

ROOT = pathlib.Path(__file__).resolve().parents[1]
import os
D = pathlib.Path(os.environ.get("MC_RES", ROOT / "out/mc_f90/res"))
Y_BON, Y_IBA = 4.930679728490592e-3, 6.9e-3
IBA_T9 = {200: 5.43, 220: 1.19, 230: 0.6}           # µSv/h contacto, KIUBE 300 doble haz 2×150 µA, concreto 2.35, frente al blanco
IBA_T8_1m = {"frente_0": 1.8e8, "90": 0.41e8}         # µSv/h a 1 m del ciclotrón (dentro de la bóveda), doble haz 2×150


def cargar(tag, cara):
    return [leer_iba.leer(f) for f in sorted(glob.glob(str(D / f"{tag}_{cara}_s*_malla.bin")))]


def interp_log(P, z, key="parche_max"):
    Z = np.array([p["prof_cm"] for p in P]); V = np.array([p[key] for p in P])
    ok = np.isfinite(V) & (V > 0); Z, V = Z[ok], V[ok]
    if z <= Z.max():
        return float(10 ** np.interp(z, Z, np.log10(V)))
    return None


def main(tag, rho):
    leer_iba.NPS = 300 * 6.2415e12 * Y_IBA                       # renormaliza a IBA
    out = {"tag": tag, "rho_corrida": rho, "rendimiento_n_p": Y_IBA, "factor_vs_BON21": Y_IBA / Y_BON, "caras": {}}
    for cara, c in A.CARAS.items():
        Ms = cargar(tag, cara)
        if not Ms:
            continue
        P = A.perfil(Ms, c); t = c["t"] / 10
        fit = A.ajustar(P, t)
        obj = A.OBJ[cara][0]
        r = {"lotes": len(Ms), "historias": float(sum(M["N"] for M in Ms)), "espesor_modelo_cm": t, "perfil": P, "ajuste": fit,
             "objetivo_uSv_h": obj}
        # requerido en la densidad de la corrida (directo, parche) y convertido a tabique 2.5 por masa
        Pin = [p for p in P if 0 < p["prof_cm"] < t and np.isfinite(p["parche_max"]) and p["parche_max"] > 0]
        zr = None
        for i, p in enumerate(Pin):
            if all(q["parche_max"] <= obj for q in Pin[i:]):
                a, b = (Pin[i - 1], p) if i else (p, p)
                if a is b:
                    zr = p["prof_cm"]
                else:
                    la, lb = math.log10(a["parche_max"]), math.log10(b["parche_max"])
                    zr = a["prof_cm"] + (math.log10(obj) - la) / (lb - la) * 10
                zr += 5
                break
        if zr is None and fit:                                     # el objetivo cae más allá del muro modelado → extrapola el ajuste
            fp = np.median([p["parche_max"] / p["H"] for p in Pin[-6:] if p["H"] > 0])
            zr = (math.log10(obj / fp) - fit["a"]) / fit["b"]; r["requerido_por_extrapolacion"] = True
        if zr is not None:
            z25 = zr * rho / 2.5
            r.update(requerido_cm_rho_corrida=round(zr, 1), requerido_cm_tabique25=round(z25, 1), modulo_13cm_tabique25=13 * math.ceil(z25 / 13))
        if fit and cara == "este":
            fp = np.median([p["parche_max"] / p["H"] for p in Pin[-8:] if p["H"] > 0])
            cmp = {}
            for zi, v in IBA_T9.items():
                zc = zi * 2.35 / rho                                  # mismo espesor másico que zi cm de concreto 2.35
                mc_dir = interp_log(Pin, zc)
                mc_fit = 10 ** (fit["a"] + fit["b"] * zc) * fp
                cmp[str(zi)] = {"IBA_uSv_h": v, "prof_equiv_cm": round(zc, 1), "MC_parche_uSv_h": mc_dir, "MC_ajuste_x_parche_uSv_h": float(mc_fit),
                                "IBA_sobre_MC_ajuste": float(v / mc_fit)}
            r["comparacion_IBA_tabla9"] = cmp
            # tasa en la bóveda a ~1 m frente al blanco T1 (x≈0.39 m, y≈0.89+1.0 m, z≈0)
            M0 = Ms[0]
            def punto(x, y, z):
                vals = []
                for M in Ms:
                    i = np.argmin(abs(M["x"] - x)); j = np.argmin(abs(M["y"] - y)); k = np.argmin(abs(M["z"] - z))
                    vals.append(float((M["Hn"] + M["Hg"])[i, j, k]))
                return float(np.mean(vals))
            r["boveda_1m_frente_T1_uSv_h"] = punto(390, 1890, 0); r["IBA_tabla8_1m"] = IBA_T8_1m
        out["caras"][cara] = r
    json.dump(out, open(ROOT / f"out/mc_f90/{tag}.json", "w"), indent=1, default=float)
    for k, r in out["caras"].items():
        print(k, r["lotes"], "TVL", r["ajuste"] and round(r["ajuste"]["TVL_cm"], 1), "req(rho)", r.get("requerido_cm_rho_corrida"),
              "req 2.5", r.get("requerido_cm_tabique25"), "mod", r.get("modulo_13cm_tabique25"))
        if "comparacion_IBA_tabla9" in r:
            for zi, v in r["comparacion_IBA_tabla9"].items():
                print("   IBA", zi, v)
            print("   1 m frente a T1:", r["boveda_1m_frente_T1_uSv_h"], "IBA", IBA_T8_1m["frente_0"])
    return out


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 2.35)
