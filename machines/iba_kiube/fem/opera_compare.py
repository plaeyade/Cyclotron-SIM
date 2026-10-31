"""Compara un mapa Bz del plano medio exportado de Opera-3d (o de cualquier código) con NGSolve y las medidas.
Uso: python fem/opera_compare.py tabla.csv [--gauss] [--mm]
La tabla debe tener columnas X, Y, (Z,) BZ (cabecera con esos nombres, separador coma o espacios)."""
import argparse, json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]


def metrics(X, Y, Bz, rs=np.arange(0.02, 0.60, 0.01)):
    r = np.hypot(X, Y); ph = np.degrees(np.arctan2(Y, X)) % 90.0
    out = []
    for r0 in rs:
        m = np.abs(r - r0) < 0.005
        if m.sum() < 20:
            continue
        b = Bz[m]; Bm = b.mean(); F = ((b - Bm) ** 2).mean() / Bm ** 2
        hill = b[np.abs(ph[m] - 45.0) < 10].mean()                      # colinas en 45° + k·90°
        vm = (ph[m] < 3) | (ph[m] > 87)                                   # valles en k·90°
        vall = b[vm].mean() if vm.any() else np.nan
        out.append((r0, Bm, F, hill, vall))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tabla"); ap.add_argument("--gauss", action="store_true")
    ap.add_argument("--mm", action="store_true", help="coordenadas en mm (por defecto m)")
    a = ap.parse_args()
    raw = open(a.tabla).read().replace(",", " ").split("\n")
    head = raw[0].upper().split(); data = np.loadtxt(raw[1:])
    col = {k: head.index(k) for k in ("X", "Y", "BZ")}
    X, Y, Bz = data[:, col["X"]], data[:, col["Y"]], data[:, col["BZ"]]
    if a.mm:
        X, Y = X / 1e3, Y / 1e3
    if a.gauss:
        Bz = Bz / 1e4
    Bz = np.abs(Bz)
    M = metrics(X, Y, Bz)
    res = {"r_m": M[:, 0].tolist(), "Bavg_T": M[:, 1].tolist(), "flutter": M[:, 2].tolist(),
           "B_colina_T": M[:, 3].tolist(), "B_valle_T": M[:, 4].tolist()}
    at = lambda k: float(np.interp(0.45, M[:, 0], M[:, k]))
    ref = {"KEY22_Bavg_T(0.38-0.47 m)": 1.3253, "KEY22_flutter": 0.394, "S09_B_colina_T": 2.15, "S09_B_valle_T": 0.45}
    ng = pathlib.Path(ROOT / "out/fem_mag/ngsolve/resumen.json")
    print(json.dumps({"tabla": {"Bavg_045": at(1), "F_045": at(2), "Bh_045": at(3), "Bv_045": at(4)},
                      "medidas": ref, "ngsolve": json.load(open(ng)) if ng.exists() else None}, indent=1))
    json.dump(res, open(pathlib.Path(a.tabla).with_suffix(".metricas.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
