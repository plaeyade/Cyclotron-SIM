"""Validación del modelo contra los documentos v0.6 (VE24, KEY22, IPAC26). Nada de esto se ajustó:
compara salidas del modelo con datos publicados. Uso: python validation/validar_documentos.py
Escribe validation/validacion_v06.json y validation/validacion_v06.png."""
import csv, json, math, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from kiube import load_params                                    # noqa: E402
from kiube.beam import extraction, variable_energy, HardEdgeField  # noqa: E402

P, _ = load_params()
E = extraction(P)
out = {}


def cross(path, R):
    r = np.hypot(path[:, 0], path[:, 1]); i = int(np.argmax(r >= R))
    return math.degrees(math.atan2(path[i, 1], path[i, 0]))


def dphi(a, hc):
    return (a - hc + 180) % 360 - 180


# 1) trayectorias 18/13 MeV vs VE24 Fig.3
ref = [r for r in csv.reader(l for l in open(ROOT / "validation/ve24_fig3_cruces.csv") if not l.startswith("#"))][1:]
V13 = variable_energy(P, E, 13.0, sel=[0, 1])
rows = []
for T, sysn, rm, d in ref:
    T, sysn, rm, d = float(T), int(sysn), float(rm), float(d)
    path = E["paths"][sysn] if T == 18 else V13["ports"][sysn]["path"]
    m = dphi(cross(path, rm * 1e3), 45.0)
    rows.append(dict(T_MeV=T, sistema=sysn, r_m=rm, ibA=d, modelo=round(m, 2), dif=round(m - d, 2)))
out["trayectorias_VE24"] = rows
out["trayectorias_VE24_max_abs_dif_deg"] = max(abs(r["dif"]) for r in rows)

# 2) energía variable: horquilla larga y cuña por energía (puerto del sistema 2)
ve = []
for T in (13.0, 14.0, 15.0, 16.0, 17.0):
    p = variable_energy(P, E, T, sel=[1])["ports"][0]
    ve.append(dict(T_MeV=T, r_strip_mm=round(p["r_strip"], 1),
                   horquilla_extra_mm=round(math.hypot(E["strippers"][1]["x"], E["strippers"][1]["y"]) - p["r_strip"], 1),
                   dphi_stripper_deg=round(p["dphi_hill_strip"] - E["strippers"][1]["dphi_hill"], 2),
                   cuna_deg=round(p["wedge_deg"], 2), cruce_820_dphi_hill=round(p["dphi_hill_820"], 2)))
out["energia_variable"] = ve
# margen dentro de la abertura del puerto a r = 820 (¿colisión con el yugo?)
Y = P["yoke"]; hw = [w / 2 for w in Y["port_widths"]]
marg = []
for sysn in (0, 1):
    p13 = variable_energy(P, E, 13.0, sel=[sysn])["ports"][0]
    c18 = math.radians(E["ports"][sysn]["phi"]); c13 = math.radians(p13["phi_820"])
    off = 820 * math.sin(c13 - c18)
    marg.append(dict(sistema=sysn, desplazamiento_mm=round(off, 1), semiancho_mm=hw[sysn],
                     cabe=bool(abs(off) + 5 < hw[sysn])))
out["energia_variable_cabe_en_puerto_13MeV"] = marg

# 3) <B> y flutter (hard-edge) vs KEY22 Fig.3
F = HardEdgeField(P)
fh = P["beam"]["hill_eff_angle"] / (360 / P["pole"]["n_sectors"])
Bh, Bv = P["beam"]["B_hill_T"], P["beam"]["B_valley_T"]
Bavg = fh * Bh + (1 - fh) * Bv
flut = fh * (1 - fh) * (Bh - Bv) ** 2 / Bavg ** 2
d = [r for r in csv.reader(l for l in open(ROOT / "validation/kiube_Bavg_flutter_KEY22.csv") if not l.startswith("#"))][1:]
rr = np.array([float(x[0]) for x in d]); BB = np.array([float(x[1]) if x[1] else np.nan for x in d])
FF = np.array([float(x[2]) if x[2] else np.nan for x in d])
sel = (rr > 0.38) & (rr < 0.47)
out["campo_KEY22"] = dict(Bavg_modelo_G=round(Bavg * 1e4, 1), Bavg_IBA_r038_047_G=round(float(np.nanmean(BB[sel])), 1),
                          flutter_modelo=round(flut, 3), flutter_IBA_r038_047=round(float(np.nanmean(FF[sel])), 3),
                          flutter_IBA_max=round(float(np.nanmax(FF)), 3))

# 4) masas y potencia de bobina vs IPAC26
CO = P["coil"]
Vc = math.pi * ((CO["inner_radius"] + CO["radial_width"]) ** 2 - CO["inner_radius"] ** 2) * CO["height"] * 1e-9
Lm = 2 * math.pi * (CO["inner_radius"] + CO["radial_width"] / 2) / 1e3
A_cu = 0.738 * CO["radial_width"] * CO["height"] * 1e-6
rho = 1.72e-8 * (1 + 0.00393 * 20)
P_coils = 2 * rho * Lm * (CO["NI_total"] / 2) ** 2 / A_cu
out["bobinas_IPAC26"] = dict(cobre_t=round(2 * Vc * 0.738 * 8.96, 2), cobre_IPAC26_t=2.0,
                             potencia_bobinas_kW_con_NI=round(P_coils / 1e3, 1), potencia_total_IPAC26_kW=19.0,
                             J_A_mm2=round(CO["NI_total"] / 2 / (A_cu * 1e6), 2))
json.dump(out, open(ROOT / "validation/validacion_v06.json", "w"), indent=1, ensure_ascii=False)

# figura
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot(rr, BB, "r-", label="IBA (KEY22 Fig.3) ⟨B⟩")
ax[0].axhline(Bavg * 1e4, color="k", ls="--", label=f"modelo hard-edge {Bavg*1e4:.0f} G")
a2 = ax[0].twinx(); a2.plot(rr, FF, "r:", label="IBA flutter"); a2.axhline(flut, color="gray", ls=":", label=f"modelo F={flut:.3f}")
a2.set_ylim(0, 0.5); a2.set_ylabel("flutter F"); ax[0].set_ylim(11500, 14000)
ax[0].set_xlabel("r [m]"); ax[0].set_ylabel("⟨Bz⟩ [G]"); ax[0].legend(loc="lower left", fontsize=8); a2.legend(loc="lower center", fontsize=8)
ax[0].set_title("Campo medio y flutter del KIUBE")
for T, c in ((18, "g"), (13, "r")):
    s = [r for r in rows if r["T_MeV"] == T]
    ax[1].plot([r["r_m"] for r in s], [r["ibA"] for r in s], c + "o", label=f"IBA {T} MeV")
    ax[1].plot([r["r_m"] for r in s], [r["modelo"] for r in s], c + "x", ms=10, label=f"modelo {T} MeV")
ax[1].set_xlabel("r [m]"); ax[1].set_ylabel("azimut desde el centro de colina [°]")
ax[1].set_title("Trayectorias de extracción (VE24 Fig.3)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(ROOT / "validation/validacion_v06.png", dpi=130)
print(json.dumps(out, indent=1, ensure_ascii=False))
