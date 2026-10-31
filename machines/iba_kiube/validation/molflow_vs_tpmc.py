"""Molflow+ 2.11.1 (CERN, molflowCLI) frente al TPMC propio (kiube/tpmc.py), en los mismos casos:
  A. Transmisión de tubos (Clausing) L/R = 1, 2, 10.
  B. Pump-down del tutorial «Time-dependent examples» (inyección 10–20 ms, 1 l/s por tapa): modo con momentos.
  C. Cámara del KIUBE (STL del LOD vacuum, 72 k facetas): H₂O de metales a 10 h y H₂ de la fuente PIG,
     cargando los MISMOS archivos --setParamsByFile que se entregan para Molflow (prueba también el flujo de trabajo).
Uso: python validation/molflow_vs_tpmc.py [--ndes 3e5] [--solo A B C]
Requiere molflowCLI (variable MOLFLOW_CLI o carpeta molflow/ del proyecto)."""
import argparse, json, math, pathlib, sys, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube.molflow import write_xml, run_cli, read_results, find_cli          # noqa: E402
from kiube.tpmc import TPMC, tube_mesh, sticking_for_speed, CLAUSING, load_gas_mesh   # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--ndes", type=float, default=3e5)
ap.add_argument("--solo", nargs="*", default=["A", "B", "C"]); a = ap.parse_args()
W = ROOT / "out/molflow_vs_tpmc"; W.mkdir(parents=True, exist_ok=True)
cli = find_cli(); print("molflowCLI:", cli)
out = json.load(open(W / "comparacion.json")) if (W / "comparacion.json").exists() else {}
out["molflow"] = "MolFlow 2.11.1 (Fedora build), molflowCLI"

# ---------------------------------------------------------------- A. Clausing
if "A" in a.solo:
    rows = []
    for LR in (1.0, 2.0, 10.0):
        m, c1, c2 = tube_mesh(R=1.0, L=LR, n=96)
        s = np.where(c1 | c2, 1.0, 0.0)
        q = np.where(c1, m.area_faces / m.area_faces[c1].sum() * 1e-6, 0.0)
        x = W / f"clausing_{LR:g}.xml"
        write_xml(x, m.triangles, sticking=s, outgassing_pa_m3_s=q, mass=28.0)
        t0 = time.time(); R = read_results(run_cli(x, ndes=a.ndes, outdir=W / f"res_clausing_{LR:g}")); tm = time.time() - t0
        M = R["moments"][0]
        Wm = M["abs"][c2].sum() / M["totalDes"]
        tp = TPMC(m, s, seed=3); r = tp.run(np.where(c1, m.area_faces, 0.0), n_particles=int(a.ndes))
        Wt = r["absorbed"][c2].sum() / r["n"]
        sig = math.sqrt(CLAUSING[LR] * (1 - CLAUSING[LR]) / a.ndes)
        rows.append(dict(L_R=LR, Clausing=CLAUSING[LR], Molflow=round(Wm, 5), TPMC=round(Wt, 5), sigma=round(sig, 5),
                         t_molflow_s=round(tm, 1)))
        print("[A]", rows[-1], flush=True)
    out["A_clausing"] = rows

# ---------------------------------------------------------------- B. pump-down del tutorial (dependiente del tiempo)
if "B" in a.solo:
    m, c1, c2 = tube_mesh(R=1.0, L=1.0, n=96)
    s = np.zeros(len(m.faces))
    for cap in (c1, c2):
        s[cap] = sticking_for_speed(1.0, m.area_faces[cap].sum(), 28.0)
    q = np.where(c1, m.area_faces / m.area_faces[c1].sum() * 1e-5, 0.0)       # 1e-5 Pa·m³/s = 1e-4 mbar·l/s
    pulse = [(0.0, 0.0), (0.0099999, 0.0), (0.01, 1.0), (0.02, 1.0), (0.020001, 0.0), (1.0, 0.0)]
    dt = 1e-4
    x = W / "tutorial_pumpdown.xml"
    write_xml(x, m.triangles, sticking=s, outgassing_pa_m3_s=q, mass=28.0, outgassing_param=np.where(c1, 0, -1),
              params={"inyeccion": pulse}, moments=[f"{dt},{dt},0.03"], time_window=dt, const_flow=False)
    Nmol = 1e-5 * 0.01 / (1.380649e-23 * 293.15)                            # moléculas reales del pulso
    t0 = time.time(); R = read_results(run_cli(x, ndes=a.ndes * 3, outdir=W / "res_tutorial"), n_molecules=Nmol,
                                       n_tp=a.ndes * 3); tm = time.time() - t0
    side = ~(c1 | c2); A = R["area_m2"]
    tt = np.arange(1, len(R["moments"])) * dt
    Pm = np.array([(M["P_mbar"][side] * A[side]).sum() / A[side].sum() for M in R["moments"][1:]])
    V = math.pi * 1e-3; tau = V / 2.0
    Pa = np.where(tt < 0.01, 0, np.where(tt < 0.02, 5e-5 * (1 - np.exp(-(tt - 0.01) / tau)),
                                         5e-5 * (1 - math.exp(-0.01 / tau)) * np.exp(-(tt - 0.02) / tau)))
    tj = json.load(open(ROOT / "validation/molflow_tutorial_pumpdown.json"))
    plat = float(Pm[(tt > 0.015) & (tt < 0.0199)].mean())
    dec = (tt > 0.0205) & (tt < 0.024) & (Pm > 0)
    tau_m = -1 / np.polyfit(tt[dec], np.log(Pm[dec]), 1)[0]
    out["B_tutorial"] = dict(P_meseta_Molflow=plat, P_meseta_TPMC=tj["P_estacionaria_mbar"], P_ref=5e-5,
                             tau_Molflow_ms=tau_m * 1e3, tau_TPMC_ms=tj["tau_TPMC_ms"], tau_analitico_ms=tau * 1e3,
                             t_molflow_s=round(tm, 1))
    np.savez(W / "tutorial_curvas.npz", t=tt, P_molflow=Pm, P_analitico=Pa)
    print("[B]", out["B_tutorial"], flush=True)

# ---------------------------------------------------------------- C. KIUBE
if "C" in a.solo:
    d = ROOT / "out/vacuum/molflow"
    m = load_gas_mesh(d / "kiube_vacio_cm.stl")
    G = json.load(open(d / "kiube_molflow_grupos.json"))["grupos"]
    T = json.load(open(ROOT / "out/vacuum/tpmc_pumpdown.json"))
    x = W / "kiube_vacio.xml"
    info = write_xml(x, m.vertices[m.faces], mass=18.015)                      # geometría sola: el resto por archivo
    c = m.triangles_center * 10.0; A_cm2 = m.area_faces
    walls = np.zeros(len(m.faces), bool)
    for g, e in G.items():
        if g.startswith("pared_"):
            walls[e["facetas"][0] - 1:e["facetas"][1]] = True
    beam = (np.hypot(c[:, 0], c[:, 1]) < 500) & (np.abs(c[:, 2]) < 70) & walls
    res = {}
    for esp, pf, mass in (("H2O_10h", "kiube_molflow_H2O_10h.txt", 18.015), ("H2_fuente", "kiube_molflow_H2_fuente.txt", 2.016)):
        t0 = time.time()
        R = read_results(run_cli(x, ndes=a.ndes, outdir=W / f"res_kiube_{esp}", params_file=d / pf, threads=2,
                                 log=W / f"log_{esp}.txt"))
        M = R["moments"][0]; P = M["P_mbar"]
        Pb = float((P * A_cm2)[beam].sum() / A_cm2[beam].sum())
        res[esp] = dict(P_haz_Molflow=Pb, t_molflow_s=round(time.time() - t0, 1), totalDes=M["totalDes"])
        print("[C]", esp, res[esp], flush=True)
    # TPMC equivalente
    R_tot = sum(v["R_haz"] for k, v in T["grupos_fuente"].items() if k.startswith("pared_"))
    res["H2O_10h"]["P_haz_TPMC"] = 3.0e-10 * R_tot
    res["H2_fuente"]["P_haz_TPMC"] = T["H2_fuente"]["P_haz_mbar"]
    for k in res:
        res[k]["Molflow_sobre_TPMC"] = res[k]["P_haz_Molflow"] / res[k]["P_haz_TPMC"]
    out["C_kiube"] = dict(facetas=info["n_facets"], zona_haz="paredes con r < 500 mm y |z| < 70 mm", **res)
json.dump(out, open(W / "comparacion.json", "w"), indent=1, ensure_ascii=False, default=float)
print(json.dumps(out, indent=1, ensure_ascii=False, default=float))
