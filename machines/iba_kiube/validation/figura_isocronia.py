"""Isocronía del imán FE (v0.33–v0.34): error de frecuencia de revolución por energía en las órbitas de equilibrio
(fem/orbitas.py) para el mapa original, las iteraciones del ajuste de colinas (fem/shim_iso.py) y el campo medido de
KEY22 (validación). Salidas: validation/isocronia.{png,json}"""
import json, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; V = ROOT / "out/fem_mag_variantes"
casos = [("KEY22 (campo medido, KIUBE real)", ROOT / "out/orbitas/equilibrio_key22.json", "#3b6fd8", 1.0038),
         ("FE original (v0.32)", V / "shim_base/orbitas.json", "#888888", 1.0),
         ("iteración 2 (perfil de colinas)", V / "shim_it2/orbitas.json", "#8e44ad", 1.0035),
         ("iteración 3 + escalón 2 a 520 mm", V / "shim_it3_L520/orbitas.json", "#16a085", 1.011),
         ("iteración 7 (no construible: invade el dee)", V / "shim_it7/orbitas.json", "#e67e22", 1.007),
         ("v0.34 adoptado: iteración 8 (con cotas)", V / "shim_it8/orbitas.json", "#d0182a", 1.0092),
         ("v0.34 FE oficial (CAD regenerado)", ROOT / "out/orbitas/equilibrio_v34.json", "#111111", 1.0098)]
ventanas = {"KEY22 (campo medido, KIUBE real)": "104° a 17 MeV · 3D ≈50°", "FE original (v0.32)": "0° (no pasa de 1.3–3.7 MeV)",
            "iteración 2 (perfil de colinas)": "90° a 17 MeV · sin órbita de 18", "iteración 3 + escalón 2 a 520 mm": "60° a 17 · 54° a 18 MeV",
            "iteración 7 (no construible: invade el dee)": "108° a 18 MeV · 3D 40°",
            "v0.34 adoptado: iteración 8 (con cotas)": "96° a 18 MeV · 3D ≈40° (80 %)",
            "v0.34 FE oficial (CAD regenerado)": "3D 36/70 (−30…+10°)"}
fig, ax = plt.subplots(figsize=(9, 5))
out = []
for lab, f, col, k in casos:
    d = json.load(open(f)); E = d["orbitas"] if isinstance(d, dict) else d
    T = [o["T_MeV"] for o in E]; e = [(o["f_rev_MHz"] * k / (40.640 / 2) - 1) * 100 for o in E]
    ax.plot(T, e, "o-", ms=3, color=col, lw=1.4, label=f"{lab} — nivel {100 * (k - 1):+.2f} % · ventana {ventanas[lab]}")
    out.append(dict(caso=lab, nivel=k, T_MeV=T, error_pct=e, nu_r=[o["nu_r"] for o in E], nu_z=[o["nu_z"] for o in E],
                    r_mm=[o["r_medio_m"] * 1e3 for o in E], ventana=ventanas[lab]))
ax.axhspan(-0.05, 0.05, color="#16a085", alpha=0.08)
ax.set_ylim(-4, 1.5); ax.set_xlabel("energía [MeV]"); ax.set_ylabel("Δf_rev / (f_RF/2) [%] con el nivel óptimo")
ax.set_title("Isocronía del imán: órbitas de equilibrio en el mapa FE (f_RF 40.640 MHz, h = 2)", fontsize=10)
ax.grid(alpha=.3); ax.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(ROOT / "validation/isocronia.png", dpi=140)
json.dump(out, open(ROOT / "validation/isocronia.json", "w"), indent=1, ensure_ascii=False)
print("ok")

# aceptancia 3D con betatrón (±6 mm radial, ±2 mm vertical): KEY22 frente al v0.34
import numpy as np
fig, ax = plt.subplots(figsize=(8.5, 4.2))
for lab, f, col in (("KEY22 (campo medido), nivel +0.38 %", ROOT / "out/orbitas/beta_key22_1.0038.json", "#3b6fd8"),
                    ("v0.34 (iteración 8), nivel +0.92 %", ROOT / "out/orbitas/beta_shim_it8_1.0092.json", "#d0182a"),
                    ("v0.34 (iteración 8), nivel +0.97 %", ROOT / "out/orbitas/beta_shim_it8_1.0097.json", "#e67e22"),
                    ("v0.34 FE oficial, nivel +0.98 % (operación)", ROOT / "out/orbitas/beta_oficial_v34_1.0098.json", "#111111")):
    if not f.exists():
        continue
    R = json.load(open(f)); fa = np.array(R["fases"]); m = np.array(R["motivo"])
    b = np.arange(-60, 11, 10); x = b[:-1] + 5
    fr = [100 * np.mean(m[(fa >= lo) & (fa < lo + 10)] == "llega") if ((fa >= lo) & (fa < lo + 10)).any() else 0 for lo in b[:-1]]
    ax.plot(x, fr, "o-", color=col, label=f"{lab}: {int((m == 'llega').sum())}/{len(m)}")
    out.append(dict(caso=lab, fases_centro=x.tolist(), llega_pct=fr))
ax.set_xlabel("fase inicial [° RF] (tras la región central, T0 = 0.75 MeV)"); ax.set_ylabel("% que llega a 18 MeV")
ax.set_title("Aceleración 3D completa hasta 18 MeV con betatrón (±6 mm radial, ±2 mm vertical)", fontsize=10)
ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_ylim(0, 105)
fig.tight_layout(); fig.savefig(ROOT / "validation/aceleracion_3d.png", dpi=140)
json.dump(out, open(ROOT / "validation/isocronia.json", "w"), indent=1, ensure_ascii=False)
