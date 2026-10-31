"""¿Hace falta SynradG4 (Geant4, Natochii et al., arXiv:2408.11709) o SynRad+ (radiación de sincrotrón) en el KIUBE? Cálculo con la física clásica (Larmor/Liénard).
Uso: python validation/evaluar_synrad.py  → validation/synrad_evaluacion.json"""
import json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube import load_params
from kiube.beam import extraction
P, _ = load_params(); E = extraction(P)
e, eps0, c, hbar_c = 1.602176634e-19, 8.8541878128e-12, 299792458.0, 197.3269804e-9   # hbar·c en eV·m
m = 939.294e6
T = P["beam"]["T_MeV"] * 1e6
g = 1 + T / m; b = math.sqrt(1 - 1 / g ** 2)
rho_h = E["rho_hill"] / 1e3
rho_v = rho_h * P["beam"]["B_hill_T"] / P["beam"]["B_valley_T"]
f_h = P["beam"]["hill_eff_angle"] / 90.0
L = 2 * math.pi * E["r_ext_hill"] / 1e3
I_ds = L * (f_h / rho_h ** 2 + (1 - f_h) / rho_v ** 2)                  # ∮ ds/ρ²
U_J = e ** 2 * b ** 3 * g ** 4 / (6 * math.pi * eps0) * I_ds            # pérdida por vuelta (Liénard)
U_eV = U_J / e
Ec_eV = 1.5 * hbar_c * g ** 3 / rho_h                                     # energía crítica
f_rev = E["f_rev_MHz"] * 1e6
I = P["beam"]["current_uA"] * 1e-6
turns = T / (P["beam"]["dE_turn_keV"] * 1e3)
N_circ = I * turns / f_rev / e
P_tot_W = N_circ * U_J * f_rev
# lo que sí importa: desorción inducida por las pérdidas de haz (H⁰ de stripping en gas) — orden de magnitud
lost = P["vacuum"]["gas_loss_fraction"] * I / e                          # iones/s que chocan con la pared
eta = (1.0, 10.0)                                                          # moléculas/ión, superficies sin hornear (EST)
Q_des = [lost * x / 2.4e19 for x in eta]                                  # mbar·l/s (1 mbar·l ≈ 2.4e19 moléculas a 20 °C)
out = dict(energia_MeV=T / 1e6, gamma=g, beta=b, rho_colina_m=rho_h, perdida_por_vuelta_eV=U_eV,
           energia_critica_eV=Ec_eV, frecuencia_fundamental_MHz=f_rev / 1e6,
           particulas_circulando=N_circ, potencia_total_radiada_W=P_tot_W,
           conclusion=("No: la radiación de un protón/H⁻ de 18 MeV (γ = 1.019) es radiación de ciclotrón de muy baja "
                       "frecuencia (armónicos de 20 MHz, energía crítica ~1e-6 eV). Pierde ~1e-10 eV por vuelta y el haz "
                       "completo radia ~1e-11 W; los fotones no pueden fotodesorber gas (umbral ~10 eV). SynRad+ es para "
                       "anillos de electrones o protones de TeV."),
           synradg4=("SynradG4 (Natochii et al., arXiv:2408.11709) extiende Geant4 para fondos de radiación de sincrotrón de "
                     "haces de ELECTRONES de GeV (ePIC/EIC): emisión con G4SynchrotronRadiation y reflexión de rayos X en las "
                     "paredes. Para H⁻/p de 18 MeV, los fotones (≈1e-6 eV, ondas de radio de 20 MHz) están 9 órdenes por "
                     "debajo del corte de producción de gammas de Geant4 (~1 keV) y no hay nada que trazar: implementarlo "
                     "daría cero. Lo que sí conviene simular en Geant4/OpenMC (ya preparado en el LOD mc): los H⁰ de "
                     "stripping en gas que golpean la pared (perdidas_h0_impactos.csv) y los neutrones de los blancos."),
           lo_que_si_importa=dict(desorcion_por_perdidas_mbar_l_s=Q_des, eta_moleculas_por_ion=eta,
                                  nota="H⁰ que chocan con la pared (15 % del haz): carga de gas 3–4 órdenes por debajo de los "
                                       "0.034 mbar·l/s de H₂ de la fuente; incluible en el TPMC como fuente en los puntos de "
                                       "perdidas_h0_impactos.csv"))
json.dump(out, open(ROOT / "validation/synrad_evaluacion.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps({k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in out.items() if k != "lo_que_si_importa"}, indent=1, ensure_ascii=False))
print("desorción:", [f"{q:.2e}" for q in Q_des], "mbar·l/s")
