"""Lectura de mallas bunker_iba (RPOS F89). H*(10) en µSv/h a 300 µA en un blanco."""
import numpy as np, glob, os


def completa(f):
    """True si la malla está entera (un lote en escritura o truncado se omite, no rompe el análisis)."""
    with open(f, "rb") as fh:
        nx, ny, nz = np.frombuffer(fh.read(12), np.int32)
    return os.path.getsize(f) == 12 + 56 + 5 * 4 * int(nx) * int(ny) * int(nz)
NPS = 300 * 6.2415e12 * 4.930679728490592e-3          # n/s
def leer(f):
    with open(f, "rb") as fh:
        nx, ny, nz = np.frombuffer(fh.read(12), np.int32)
        g = np.frombuffer(fh.read(56), np.float64)
        n = nx * ny * nz
        A = []
        for _ in range(5):
            b = fh.read(4 * n)
            if len(b) != 4 * n:
                raise ValueError(f"malla incompleta: {f}")
            A.append(np.frombuffer(b, np.float32).reshape(nx, ny, nz).astype(np.float64))
    N = g[4]; vol = (g[3] / 10) ** 3
    k = NPS * 3600 * 1e-6 / vol / N
    x = g[0] + g[3] * (np.arange(nx) + 0.5); y = g[1] + g[3] * (np.arange(ny) + 0.5); z = g[2] + g[3] * (np.arange(nz) + 0.5)
    return dict(x=x, y=y, z=z, Hn=A[0] * k, Hg=A[1] * k, N=N)
