#!/usr/bin/env python3
"""OpenMC (≥0.14, compilado con DAGMC) — neutrones de los blancos del Cyclone KIUBE.

  python run_openmc.py --targets T1 T2 --particles 200000 --batches 20

Resultados por segundo de operación (normalizados a la intensidad de los blancos activos):
  - flujo de neutrones en malla (térmico / epitérmico / rápido) y tasa de dosis H*(10) (ICRP-74/116)
  - tasas de reacción de activación por material: 63Cu(n,α)60Co y 63Cu(n,γ)64Cu en bobinas,
    109Ag(n,γ)110mAg, 59Co(n,γ)60Co en acero, 151Eu(n,γ)152Eu en concreto
Unidades: la geometría DAGMC está en cm.
"""
import argparse, json, pathlib
import numpy as np
import yaml
import openmc

HERE = pathlib.Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--targets", nargs="+", default=["T1", "T2"], help="blancos activos (KIUBE: 2 simultáneos)")
ap.add_argument("--particles", type=int, default=200000)
ap.add_argument("--batches", type=int, default=20)
ap.add_argument("--mesh-cm", type=float, default=5.0)
ap.add_argument("--bonvin-correction", action="store_true",
                help="multiplica la intensidad por 2.5 (corrección empírica de Bonvin 2021 sobre FLUKA)")
a = ap.parse_args()

M = yaml.safe_load(open(HERE / "materials.yaml"))
S = json.load(open(HERE / "sources.json"))

# ---------------- materiales (nombres = tags DAGMC)
mats = openmc.Materials()
for name, m in M.items():
    if not isinstance(m, dict) or "density" not in m or name in ("vacuum",):
        continue
    mat = openmc.Material(name=name)
    if "isotopic" in m:
        iso = m["isotopic"]
        mat.add_element("H", iso.get("H1", 2), "ao")
        if iso.get("O16"):
            mat.add_nuclide("O16", iso["O16"], "ao")
        if iso.get("O18"):
            mat.add_nuclide("O18", iso["O18"], "ao")
    else:
        enr = m.get("enriquecimiento", {})
        for el, w in m["wt"].items():
            if w <= 0:
                continue
            if el in enr:                       # isótopos en fracción atómica dentro del elemento → fracción másica
                ab = enr[el]; A = {n: openmc.data.atomic_mass(n) for n in ab}
                Mel = sum(ab[n] * A[n] for n in ab)
                for n in ab:
                    mat.add_nuclide(n, w / 100.0 * ab[n] * A[n] / Mel, "wo")
            else:
                mat.add_element(el, w / 100.0, "wo")
    mat.set_density("g/cm3", m["density"])
    mats.append(mat)
mats.export_to_xml()

# ---------------- geometría DAGMC (aire = implicit complement; 'vacuum' = vacío)
dag = openmc.DAGMCUniverse(str(HERE / "kiube_openmc.h5m"), auto_geom_ids=True)
geometry = openmc.Geometry(dag.bounded_universe(padding_distance=10.0))
geometry.export_to_xml()

# ---------------- fuente: espectro H2-18O (Bonvin 2021 Fig.7), isótropa en la cavidad
sp = S["spectrum"]
E = np.array(sp["E_MeV"]) * 1e6
EdNdE = np.array(sp["EdNdE_per_proton"])
o = np.argsort(E); E, EdNdE = E[o], EdNdE[o]
pdf = EdNdE / E                                   # dN/dE
srcs = []
k = 2.5 if a.bonvin_correction else 1.0
for s in S["sources"]:
    if s["name"] not in a.targets:
        continue
    src = openmc.IndependentSource()
    src.space = openmc.stats.Point([c / 10.0 for c in s["pos_mm"]])
    src.angle = openmc.stats.Isotropic()
    src.energy = openmc.stats.Tabular(E, pdf, interpolation="linear-linear")
    src.strength = k * (s["neutrons_per_s_if_active"] or 0.0)
    srcs.append(src)
total_strength = sum(s.strength for s in srcs)

settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = srcs
settings.particles = a.particles
settings.batches = a.batches
settings.photon_transport = False
settings.export_to_xml()

# ---------------- tallies
tallies = openmc.Tallies()
mesh = openmc.RegularMesh()
mesh.lower_left = [-450.0, -450.0, -150.0]
mesh.upper_right = [450.0, 450.0, 150.0]
mesh.dimension = [int(900 / a.mesh_cm), int(900 / a.mesh_cm), int(300 / a.mesh_cm)]
mf = openmc.MeshFilter(mesh)
ef = openmc.EnergyFilter([1e-5, 0.625, 1e5, 20e6])          # térmico / epitérmico / rápido
t = openmc.Tally(name="flux_mesh"); t.filters = [mf, ef]; t.scores = ["flux"]; tallies.append(t)
e_d, c_d = openmc.data.dose_coefficients("neutron", "AP")     # pSv·cm²
t = openmc.Tally(name="dose_mesh"); t.filters = [mf, openmc.EnergyFunctionFilter(e_d, c_d)]; t.scores = ["flux"]
tallies.append(t)
by_name = {m.name: m for m in mats}
act = [("copper_coil_homog", "Cu63", "(n,a)"), ("copper_coil_homog", "Cu63", "(n,gamma)"),
       ("copper_coil_homog", "Ag109", "(n,gamma)"), ("steel_1010", "Co59", "(n,gamma)"),
       ("concrete_portland", "Eu151", "(n,gamma)"), ("concrete_portland", "Co59", "(n,gamma)")]
for mname, nuc, rx in act:
    t = openmc.Tally(name=f"{mname}:{nuc}{rx}")
    t.filters = [openmc.MaterialFilter([by_name[mname]])]
    t.nuclides = [nuc]
    t.scores = [rx]
    tallies.append(t)
tallies.export_to_xml()

sp_file = openmc.run()
with openmc.StatePoint(sp_file) as spt:
    out = {"source_neutrons_per_s": total_strength, "targets": a.targets, "reaction_rates_per_s": {}}
    for mname, nuc, rx in act:
        tt = spt.get_tally(name=f"{mname}:{nuc}{rx}")
        out["reaction_rates_per_s"][f"{mname}:{nuc}{rx}"] = {
            "mean": float(tt.mean.sum() * total_strength), "std": float(np.sqrt((tt.std_dev ** 2).sum()) * total_strength)}
    vox = a.mesh_cm ** 3
    d = spt.get_tally(name="dose_mesh").mean.reshape(mesh.dimension[::-1])
    out["max_dose_rate_Sv_per_h_in_mesh"] = float(d.max() / vox * total_strength * 1e-12 * 3600)
json.dump(out, open("kiube_openmc_results.json", "w"), indent=1)
print(json.dumps(out, indent=1))
