"""Flujo térmico medio en el aire de la bóveda (F90) → contrato `cyclotron-sim/flujo_aire@1`.

Insumo de RadShieldMx (activación de Ar-41 en el aire, NCRP 144). Promedia el flujo de neutrones
térmicos (< 0,5 eV; longitud de traza de la malla bunker_iba) sobre el volumen interior de la
bóveda. Excluye la caja de la máquina, donde los vóxeles son hierro y no aire.
Normalización: corriente total y rendimiento IBA 6,9E-3 n/p (SRG §2.1.8.4).
Uso: MC_RES=<dir> python flujo_aire_f90.py <tag> [uA] [salida.json]
"""
import glob, json, os, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = pathlib.Path(os.environ.get("MC_RES", ROOT / "out/mc_f90/res"))
Y_IBA = 6.9e-3
INTERIOR = {"x": (-1970, 1970), "y": (-2000, 2000), "z": (-1120, 1380)}   # mm, caras internas
MAQUINA = {"x": (-1100, 1100), "y": (-1100, 1100), "z": (-1120, 1300)}    # caja conservadora del KIUBE


def leer(f):
    with open(f, "rb") as fh:
        nx, ny, nz = np.frombuffer(fh.read(12), np.int32)
        g = np.frombuffer(fh.read(56), np.float64)
        n = int(nx) * int(ny) * int(nz)
        A = []
        for _ in range(3):
            b = fh.read(4 * n)
            if len(b) != 4 * n:
                raise ValueError(f"malla incompleta: {f}")
            A.append(np.frombuffer(b, np.float32).reshape(nx, ny, nz))
    vol = (g[3] / 10) ** 3
    ejes = [g[i] + g[3] * (np.arange(m) + 0.5) for i, m in enumerate((nx, ny, nz))]
    return ejes, A[2].astype(np.float64) / 10 / vol / g[4]      # fluencia térmica por neutrón fuente [1/cm²]


def mascara(ejes, caja, dentro=True):
    m = np.ones([len(e) for e in ejes], bool)
    for i, k in enumerate("xyz"):
        sel = (ejes[i] > caja[k][0]) & (ejes[i] < caja[k][1])
        sh = [1, 1, 1]; sh[i] = len(ejes[i])
        m &= sel.reshape(sh)
    return m if dentro else ~m


def main(tag, uA=300.0, salida=None):
    S = uA * 1e-6 / 1.602176634e-19 * Y_IBA
    vals, lotes = [], 0
    for f in sorted(glob.glob(str(D / f"{tag}_*_s*_malla.bin"))):
        try:
            ejes, phi = leer(f)
        except ValueError:
            continue
        aire = mascara(ejes, INTERIOR) & mascara(ejes, MAQUINA, dentro=False)
        vals.append(phi[aire] * S)
        lotes += 1
    v = np.mean(np.stack(vals), axis=0)
    out = {
        "contrato": "cyclotron-sim/flujo_aire@1",
        "maquina": "IBA Cyclone KIUBE 300", "configuracion": tag, "corriente_uA": uA,
        "rendimiento_n_p": Y_IBA, "lotes": lotes,
        "phi_termico_aire_n_cm2_s": {"media": float(v.mean()), "mediana": float(np.median(v)),
                                       "p95": float(np.percentile(v, 95)), "max": float(v.max())},
        "volumen_promediado": {"interior_mm": INTERIOR, "excluida_maquina_mm": MAQUINA},
        "uso": "RadShieldMx ventilación Ar-41: usar la media; p95 como cota conservadora",
        "fuente": "Geant4 11.4.3 QGSP_BIC_HP, bunker_iba F90 (sin puertas de blanco)",
    }
    p = pathlib.Path(salida or ROOT / "export/flujo_aire_boveda.json")
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(out["phi_termico_aire_n_cm2_s"]), lotes)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], float(a[1]) if len(a) > 1 else 300.0, a[2] if len(a) > 2 else None)
