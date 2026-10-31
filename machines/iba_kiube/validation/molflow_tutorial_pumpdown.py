"""Validación del TPMC en modo dependiente del tiempo contra el caso de pump-down del tutorial
oficial de Molflow+ («Time-dependent examples», Kersevan & Ady, CERN, actualizado 2021 para 2.8.5):
  tubo L/R = 1 (Test > Pipe), inyección de 1e-4 mbar·l/s desde la tapa 1 entre t = 10 y 20 ms,
  1 l/s de bombeo en cada tapa → presión estacionaria 5e-5 mbar; subida y bajada exponenciales.
Aquí R = 1 cm (supuesto; el tutorial no da el radio, y τ = V/S escala con él). Gas N₂ a 20 °C.
Uso: python validation/molflow_tutorial_pumpdown.py"""
import json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube.tpmc import TPMC, tube_mesh, sticking_for_speed, vbar   # noqa: E402

R, L, M, T = 1.0, 1.0, 28.0, 293.15
mesh, cap1, cap2 = tube_mesh(R=R, L=L, n=96)
s = np.zeros(len(mesh.faces))
for cap in (cap1, cap2):
    s[cap] = sticking_for_speed(1.0, mesh.area_faces[cap].sum(), M, T)
src = cap1.astype(float) * mesh.area_faces
det = ~(cap1 | cap2)                                       # pared lateral = «vacuómetro»
dt = 0.05e-3
edges = np.arange(0, 30e-3 + dt / 2, dt)
tp = TPMC(mesh, s, seed=5)
N = 200000
res = tp.run_tof(src, det, edges, M, T, n_particles=N, batch=50000, length_scale_m=1e-2)
h = res["H"] / N                                            # impactos en el detector por molécula inyectada
A_det = mesh.area_faces[det].sum() * 1e-4                   # m²
kT = 1.380649e-23 * T; m = M * 1.66053907e-27
t = 0.5 * (edges[1:] + edges[:-1])
Q = np.where((t >= 10e-3) & (t < 20e-3), 1e-4, 0.0)        # mbar·l/s
Ndot = Q * 0.1 / kT                                         # moléculas/s  (1 mbar·l = 0.1 Pa·m³)
Gam = np.convolve(Ndot, h)[:len(t)] / A_det                 # impactos /(m² s) en el detector
P = Gam * math.sqrt(2 * math.pi * m * kT) / 100.0           # Pa → mbar
# analítico (volumen bien mezclado)
V = math.pi * R ** 2 * L * 1e-3                             # litros
S = 2.0                                                     # l/s
tau = V / S
Pa = np.where(t < 10e-3, 0.0, np.where(t < 20e-3, 5e-5 * (1 - np.exp(-(t - 10e-3) / tau)),
                                        5e-5 * (1 - math.exp(-10e-3 / tau)) * np.exp(-(t - 20e-3) / tau)))
plateau = P[(t > 15e-3) & (t < 20e-3)].mean()
# constante de tiempo efectiva medida en la bajada
dec = (t > 20.5e-3) & (t < 24e-3) & (P > 0)
tau_mc = -1 / np.polyfit(t[dec], np.log(P[dec]), 1)[0]
out = dict(caso="Molflow+ tutorial 'Time-dependent examples', ejemplo 1 (tubo L/R=1, inyección 10–20 ms)",
           R_cm=R, gas="N2", P_estacionaria_mbar=plateau, P_referencia_mbar=5e-5, dif_rel=plateau / 5e-5 - 1,
           tau_analitico_ms=tau * 1e3, tau_TPMC_ms=tau_mc * 1e3, particulas=N, vbar_m_s=vbar(M, T))
json.dump(out, open(ROOT / "validation/molflow_tutorial_pumpdown.json", "w"), indent=1, ensure_ascii=False)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(t * 1e3, P, color="#1f6fb2", lw=1.6, label="TPMC (tiempo de vuelo)")
ax.plot(t * 1e3, Pa, "k--", lw=1, label=f"analítico V/S, τ = {tau*1e3:.2f} ms")
ax.axhline(5e-5, color="#c0392b", lw=0.8, ls=":", label="estacionario del tutorial: 5·10⁻⁵ mbar")
ax.set_xlabel("t [ms]"); ax.set_ylabel("presión en la pared [mbar]"); ax.set_xlim(0, 30)
ax.set_title("Validación: caso pump-down del tutorial de Molflow+ (inyección 10–20 ms)")
ax.legend(fontsize=8); ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(ROOT / "validation/molflow_tutorial_pumpdown.png", dpi=130)
print(json.dumps(out, indent=1, ensure_ascii=False))
