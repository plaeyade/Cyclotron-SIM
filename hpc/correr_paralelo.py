#!/usr/bin/env python3
"""Lanza lotes independientes de bunker_iba en todos los núcleos (cada lote = 1 proceso, semilla distinta).
Uso: correr_paralelo.py <cara|todas> <n_nucleos (0 = todos)> [historias_por_lote=200] [lotes_por_cara=0 → 2×núcleos]
     [config=dual235.txt] [gdml=gdml/kiube_sinpuertas.gdml]
Las salidas (*_malla.bin, *_regiones.txt) van a /salida/<config>_<cara>_s<semilla>_*. Se pueden interrumpir y reanudar:
un lote con regiones.txt ya existente se salta."""
import os, subprocess, sys, pathlib, itertools, concurrent.futures as cf, time

BASE = pathlib.Path(__file__).resolve().parent
EXE = BASE / "mc_benchmark/bunker_iba"
MC = BASE / "mc_f90"
SAL = pathlib.Path(os.environ.get("SALIDA", "/salida"))
CARAS = {"este": "3.0 1 1 2000 25", "oeste": "3.0 1 -1 2000 27", "sur": "3.0 0 1 1970 28", "techo": "3.0 2 1 1380 26"}


def lote(cfg, gdml, cara, n, seed):
    tag = pathlib.Path(cfg).stem
    pre = SAL / f"{tag}_{cara}_s{seed}"
    if (SAL / f"{pre.name}_regiones.txt").exists():
        return pre.name, "ya"
    t0 = time.time()
    cmd = [str(EXE), str(MC / cfg), str(MC / gdml), str(n), str(seed), str(pre) + "_tmp", *CARAS[cara].split()]
    with open(str(pre) + ".log", "w") as log:
        r = subprocess.run(cmd, cwd=MC, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode == 0:
        os.replace(str(pre) + "_tmp_malla.bin", str(pre) + "_malla.bin")
        os.replace(str(pre) + "_tmp_regiones.txt", str(pre) + "_regiones.txt")
        os.remove(str(pre) + ".log")
    return pre.name, f"rc={r.returncode} {time.time() - t0:.0f}s"


def main():
    cara = sys.argv[1] if len(sys.argv) > 1 else "todas"
    nuc = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    nuc = nuc or os.cpu_count()
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    lotes = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    cfg = sys.argv[5] if len(sys.argv) > 5 else "dual235.txt"
    gdml = sys.argv[6] if len(sys.argv) > 6 else "gdml/kiube_sinpuertas.gdml"
    caras = list(CARAS) if cara == "todas" else [cara]
    lotes = lotes or 2 * nuc
    SAL.mkdir(parents=True, exist_ok=True)
    trabajos = [(cfg, gdml, c, n, 1000 + 100 * i + k) for i, c in enumerate(caras) for k in range(lotes)]
    trabajos.sort(key=lambda t: t[4] % 100)          # intercala caras para que todas avancen a la vez
    print(f"{len(trabajos)} lotes en {nuc} núcleos → {SAL}", flush=True)
    with cf.ProcessPoolExecutor(max_workers=nuc) as ex:
        futuros = [ex.submit(lote, *t) for t in trabajos]
        for f in cf.as_completed(futuros):
            nom, st = f.result()
            print(nom, st, flush=True)


if __name__ == "__main__":
    main()
