"""Activación del concreto de la bóveda (F90): actividad específica al final de la vida útil vs profundidad, por cara,
y espesor de la capa interior de desmantelamiento (donde Σ A_i / NL_i ≤ 1, niveles de dispensa IAEA GSR Part 3 tabla I.2).
Fluencias por grupo de la malla de bunker_iba (térmico < 0.5 eV, epitérmico 0.5 eV–0.1 MeV, rápido); reacciones (n,γ).
Uso: MC_RES=<dir> python activacion_f90.py <tag> [rho] [anios] [h_por_anio] [uA]"""
import glob, json, math, os, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = pathlib.Path(os.environ.get("MC_RES", ROOT / "out/mc_f90/res"))
Y_IBA = 6.9e-3
NA = 6.02214076e23
B = 1e-24
LN_EPI = math.log(0.1e6 / 0.5)          # letargía del grupo epitérmico
SIG_TH = math.sqrt(math.pi) / 2          # Maxwelliana a 293 K vs 2200 m/s
# (elemento padre, fracción másica en el concreto, A, abundancia isotópica, σ0 [b], RI [b], T½ [a], NL [Bq/g])
NUC = {
    "Eu-152": ("Eu", 1e-6, 151.96, 0.478, 5900.0, 3300.0, 13.54, 0.1),
    "Eu-154": ("Eu", 1e-6, 151.96, 0.522, 312.0, 1600.0, 8.59, 0.1),
    "Co-60": ("Co", 1e-5, 58.93, 1.0, 37.2, 75.0, 5.27, 0.1),
    "Ba-133": ("Ba", 0.0, 137.33, 0.00101, 6.5, 24.0, 10.55, 0.1),
    "Fe-55": ("Fe", 0.014, 55.85, 0.0585, 2.25, 1.2, 2.74, 1000.0),
}
CARAS = {"este": (1, 1, 2000), "oeste": (1, -1, 2000), "sur": (0, 1, 1970), "techo": (2, 1, 1380)}
PERP = {"este": {"x": (-1970, 1970), "z": (-1120, 1380)}, "oeste": {"x": (-1970, 1970), "z": (-1120, 1380)},
        "sur": {"y": (-2000, 2000), "z": (-1120, 1380)}, "techo": {"x": (-1970, 1970), "y": (-2000, 2000)}}


def leer(f):
    with open(f, "rb") as fh:
        nx, ny, nz = np.frombuffer(fh.read(12), np.int32)
        g = np.frombuffer(fh.read(56), np.float64)
        n = nx * ny * nz
        A = [np.frombuffer(fh.read(4 * n), np.float32).reshape(nx, ny, nz).astype(np.float64) for _ in range(5)]
    N = g[4]; vol = (g[3] / 10) ** 3
    x = g[0] + g[3] * (np.arange(nx) + 0.5); y = g[1] + g[3] * (np.arange(ny) + 0.5); z = g[2] + g[3] * (np.arange(nz) + 0.5)
    # fluencia por neutrón fuente [1/cm²]: longitud de traza (mm → cm) / volumen del vóxel
    return dict(x=x, y=y, z=z, th=A[2] / 10 / vol / N, ep=A[3] / 10 / vol / N, fa=A[4] / 10 / vol / N)


def main(tag, rho=2.35, anios=20.0, h=2000.0, uA=300.0, ba=0.0):
    NUC["Ba-133"] = NUC["Ba-133"][:1] + (ba,) + NUC["Ba-133"][2:]
    S = uA * 1e-6 / 1.602176634e-19 * Y_IBA * h / 8766.0           # n/s promedio en el año
    out = {"tag": tag, "rho": rho, "anios": anios, "h_por_anio": h, "uA": uA, "Ba_frac": ba, "caras": {}}
    for cara, (ax, sg, p0) in CARAS.items():
        fs = sorted(glob.glob(str(D / f"{tag}_{cara}_s*_malla.bin")))
        if not fs:
            continue
        Ms = [leer(f) for f in fs]; M0 = Ms[0]; c = {"x": M0["x"], "y": M0["y"], "z": M0["z"]}
        sel = [np.ones_like(c[k], bool) for k in "xyz"]
        for k, (a, b) in PERP[cara].items():
            i = "xyz".index(k); sel[i] = (c[k] > a) & (c[k] < b)
        prof = sg * c["xyz"[ax]] - p0
        filas = []
        for kk in range(0, 15):
            m = (prof > kk * 100) & (prof < (kk + 1) * 100)
            s = list(sel); s[ax] = s[ax] & m
            g = {k: np.mean([M[k][np.ix_(*s)].mean() for M in Ms]) * S for k in ("th", "ep", "fa")}   # n/cm²/s promedio
            gmax = {k: np.mean([M[k][np.ix_(*s)].max() for M in Ms]) * S for k in ("th", "ep")}
            fila = {"prof_cm": kk * 10 + 5, "phi_th": g["th"], "phi_ep": g["ep"], "phi_fa": g["fa"], "act_Bq_g": {}, "act_Bq_g_max": {}}
            suma = suma_max = 0.0
            for nom, (el, w, A, ab, s0, ri, t12, nl) in NUC.items():
                Nper_g = w * NA / A * ab                                     # átomos blanco por g de concreto
                for key, gg in (("act_Bq_g", g), ("act_Bq_g_max", gmax)):
                    R = Nper_g * (s0 * SIG_TH * gg["th"] + ri * gg["ep"] / LN_EPI) * B      # reacciones /s /g
                    lam = math.log(2) / t12
                    a = R * (1 - math.exp(-lam * anios))
                    fila[key][nom] = a
                suma += fila["act_Bq_g"][nom] / nl; suma_max += fila["act_Bq_g_max"][nom] / nl
            fila["suma_sobre_NL"] = suma; fila["suma_sobre_NL_max"] = suma_max
            filas.append(fila)
        capa = next((f["prof_cm"] + 5 for f in filas if f["suma_sobre_NL_max"] <= 1), None)
        out["caras"][cara] = {"lotes": len(Ms), "perfil": filas, "capa_desmantelamiento_cm": capa,
                              "capa_desmantelamiento_cm_tabique25": round(capa * rho / 2.5) if capa else None}
        print(cara, len(Ms), "capa", capa, "cm;", " ".join(f"{f['prof_cm']}:{f['suma_sobre_NL_max']:.3g}" for f in filas[:10]))
    OUT = pathlib.Path(os.environ.get("MC_OUT", ROOT / "out/mc_f90")); OUT.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(OUT / f"activacion_{tag}_Ba{ba}.json", "w"), indent=1, default=float)
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], *(float(v) for v in a[1:]))
