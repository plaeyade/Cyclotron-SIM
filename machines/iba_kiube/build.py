#!/usr/bin/env python3
"""Construye el modelo KIUBE para un nivel de detalle (LOD) y exporta.

  python build.py --lod neutronics   # todo (láminas µm, He, H2-18O) → STEP por material + script Cubit/DAGMC
  python build.py --lod fem_struct   # hardware ≥ 0.5 mm, un STEP por pieza → Abaqus / ANSYS Mechanical / COMSOL SM
  python build.py --lod fem_mag      # hierro + bobinas + aire conforme, mitad sup., sector 90° y 360° → COMSOL AC/DC, Opera, Maxwell
  python build.py --lod print        # piezas ≥ 0.5 mm → STL por pieza
  python build.py --lod print3d --escala 0.1 --impresora fdm   # cuerpos imprimibles → Blender (blender/)
  python build.py --lod mc           # DAGMC (OpenMC, DAG-MCNP) + GDML (Geant4) desde una malla conforme
"""
import argparse, json, sys, itertools, pathlib, math
import numpy as np
import yaml
import cadquery as cq

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from kiube import load_params, ROOT
from kiube.beam import extraction
from kiube.geometry import build, build_fem_mag

COLORS = {"steel_1010": (0.55, 0.57, 0.6), "copper_coil_homog": (0.85, 0.6, 0.2),
          "copper_ofe": (0.8, 0.35, 0.2), "copper_dee_caja": (0.8, 0.35, 0.2), "stainless_304L": (0.75, 0.75, 0.78),
          "aluminum_6061": (0.85, 0.85, 0.9), "niobium": (0.5, 0.5, 0.65),
          "havar": (0.3, 0.3, 0.35), "titanium": (0.6, 0.6, 0.7), "helium": (0.9, 0.9, 1.0),
          "water_o18": (0.2, 0.4, 0.9), "water_nat": (0.3, 0.6, 1.0), "n2_20bar": (0.4, 0.8, 0.8),
          "graphite_foil": (0.1, 0.1, 0.1), "air": (0.8, 0.9, 1.0),
          "brass_cuzn39": (0.82, 0.73, 0.51), "polymer_red": (0.85, 0.2, 0.15), "hose_pu": (0.25, 0.59, 0.89),
          "hose_epdm": (0.05, 0.05, 0.06), "cabinet_equiv": (0.92, 0.92, 0.92), "transmitter_equiv": (0.25, 0.52, 0.78),
          "polymer_black": (0.16, 0.16, 0.18), "steel_dark": (0.22, 0.23, 0.24), "paint_white": (0.92, 0.92, 0.92),
          "hdpe": (0.08, 0.08, 0.09)}


def overlaps(parts, tol=1e-3):
    bad, bbs = [], [p["shape"].BoundingBox() for p in parts]
    for i, j in itertools.combinations(range(len(parts)), 2):
        a, b = bbs[i], bbs[j]
        if a.xmax < b.xmin or b.xmax < a.xmin or a.ymax < b.ymin or b.ymax < a.ymin \
                or a.zmax < b.zmin or b.zmax < a.zmin:
            continue
        v = parts[i]["shape"].intersect(parts[j]["shape"]).Volume()
        if v > tol:
            bad.append((parts[i]["name"], parts[j]["name"], round(v, 4)))
    return bad


def save_assembly(parts, path, name):
    asm = cq.Assembly(name=name)
    for p in parts:
        asm.add(p["shape"], name=p["name"], color=cq.Color(*COLORS.get(p["mat"], (0.7, 0.7, 0.7))))
    asm.export(str(path)) if hasattr(asm, "export") else asm.save(str(path))


def mass_report(parts, mats):
    by_mat, fam = {}, {}
    for p in parts:
        m = p["shape"].Volume() * 1e-6 * mats[p["mat"]]["density"]  # mm3 → kg (g/cm3·1e-3)
        p["mass"] = m
        by_mat[p["mat"]] = by_mat.get(p["mat"], 0) + m
        f = p["name"].split("_")[0]
        fam[f] = fam.get(f, 0) + m
    return {k: round(v, 2) for k, v in by_mat.items()}, {k: round(v, 1) for k, v in fam.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lod", default="neutronics", choices=["print", "fem_struct", "fem_mag", "fem_rf", "neutronics", "mc", "print3d", "vacuum", "render", "stl_hd"])
    ap.add_argument("--escala", type=float, default=0.1, help="print3d: 0.1 = 1:10")
    ap.add_argument("--impresora", default="fdm", choices=["fdm", "resina"])
    ap.add_argument("--calidad", default="alta", choices=["media", "alta", "ultra"],
                    help="stl_hd: flecha 0.02 / 0.005 / 0.002 mm y ángulo 0.10 / 0.05 / 0.02 rad (a la escala de salida)")
    ap.add_argument("--remallar", type=float, default=None, help="stl_hd: remallado isotrópico con arista objetivo [mm]")
    ap.add_argument("--paleta", default="foto", choices=["foto", "render", "didactico"],
                    help="colores (palette.yaml, medidos en las fichas): foto = máquina real, render = render oficial, "
                         "didactico = por material físico (v0.9)")
    ap.add_argument("--sin_completo", action="store_true", help="stl_hd: no generar la unión completa")
    ap.add_argument("--puente", type=float, default=1.2, help="stl_hd: radio [mm] de los puntales que unen piezas sueltas al bloque (0 = sin puentes)")
    ap.add_argument("--molflow", default=None, help="vacuum: CSV de MolflowCLI (--writeFacetQuantities) → P(r) real; 'auto' = correr molflowCLI aquí")
    ap.add_argument("--tpmc", type=int, default=30000, help="vacuum: partículas TPMC por fuente (0 = sin TPMC)")
    ap.add_argument("--out", default=str(ROOT / "out"))
    ap.add_argument("--no-check", action="store_true")
    a = ap.parse_args()
    P, raw = load_params()
    mats = yaml.safe_load(open(ROOT / "materials.yaml"))
    E = extraction(P)
    out = pathlib.Path(a.out) / a.lod
    out.mkdir(parents=True, exist_ok=True)
    optics = {k: round(float(E[k]), 3) for k in ("hill_eff_angle", "f_rev_MHz", "r_ext_hill", "r_ext_valley",
                                                 "rho_hill", "edge_d0", "edge_R", "edge_L", "Rh_over_Lh")}
    report = {"lod": a.lod, "optics": optics,
              "strippers_dphi_hill": [round(s["dphi_hill"], 2) for s in E["strippers"][:2]],
              "edge_crossings_dphi_hill": [round(e["dphi_hill"], 2) for e in E["edge_crossings"][:2]]}

    if a.lod == "vacuum":
        from kiube.mc import mc_volumes
        from kiube.vacuum import classify_vacuum_surfaces, export_molflow, analytic_pressure, hminus_losses, lorentz_lifetime_s
        vols, _ = mc_volumes(P, E, "openmc")
        groups, meta = classify_vacuum_surfaces(vols, P)
        info = export_molflow(P, groups, meta, out / "molflow")
        pa = analytic_pressure(P, info)
        V_ = P["vacuum"]
        if a.tpmc > 0:
            from kiube.pumpdown import run as run_tpmc, plot as plot_pd
            T = run_tpmc(P, out / "molflow", n_particles=a.tpmc)
            json.dump(T, open(out / "tpmc_pumpdown.json", "w"), indent=1)
            plot_pd(P, T, out / "pumpdown.png")
            report["tpmc"] = {k: v for k, v in T.items() if k not in ("pumpdown", "H2_perfil_radial")}
        if a.molflow == "auto":                           # corre molflowCLI aquí mismo (H₂ de la fuente)
            from kiube.molflow import run_cli, find_cli
            if find_cli():
                res = run_cli(out / "molflow" / "kiube_vacio.zip", ndes=3e5, outdir=out / "molflow_resultados",
                              params_file=out / "molflow" / "kiube_molflow_H2_fuente.txt",
                              log=out / "molflow_resultados.log")
                a.molflow = str(res.parent / "facet_physics.csv")
                res.unlink()                              # XML de resultados (~100 MB): queda el .zip y el CSV
                report["molflow_cli"] = dict(resultado=str(res.relative_to(out)), particulas=3e5)
            else:
                print("[aviso] molflowCLI no encontrado: se usa el TPMC", flush=True); a.molflow = None
        if a.molflow:
            from kiube.vacuum import pressure_profile_from_molflow
            Pf, prof = pressure_profile_from_molflow(a.molflow, out / "molflow" / "kiube_molflow_grupos.json",
                                                     out / "molflow" / "kiube_vacio_cm.stl")
            pa["fuente_presion"] = "Molflow+ (perfil radial cerca del plano medio)"
            pa["perfil_molflow"] = prof
            L = hminus_losses(P, E, Pf)
        elif a.tpmc > 0:
            prof = np.array([p for p in T["H2_perfil_radial"] if p[1] == p[1]])
            Pbase = T["comparacion"]["P_24h"]
            pa["fuente_presion"] = "TPMC: perfil radial de H₂ de la fuente + H₂O a 24 h (forma; magnitud calibrada a la transmisión)"
            L = hminus_losses(P, E, lambda r_mm: float(np.interp(r_mm, prof[:, 0], prof[:, 1])) + Pbase)
        else:
            pa["fuente_presion"] = "estimador analítico Q/S"
            L = hminus_losses(P, E, pa["P_H2_operacion_mbar"] + pa["P_H2O_10h_mbar"])
        tau, Ef = lorentz_lifetime_s(P["beam"]["T_MeV"], P["beam"]["B_hill_T"])
        import shutil
        shutil.copy(ROOT / "vacuum_templates" / "README_MOLFLOW.md", out / "README_MOLFLOW.md")
        hits = L.pop("impactos")
        np.savetxt(out / "perdidas_h0_impactos.csv", hits, delimiter=",",
                   header="E_MeV,x_mm,y_mm,cos_incidencia  (H0 → pared interior de la cámara, z = 0)", comments="")
        report.update(grupos={g: v["area_cm2"] for g, v in meta.items()}, presion=pa, perdidas_gas=L,
                      lorentz={"campo_propio_V_m": Ef, "vida_propia_s": tau,
                               "fraccion_perdida_18MeV_por_m": 1 - math.exp(-1.0 / (2.99792458e8 * 0.194 * 1.019 * tau))})
    elif a.lod == "print3d":
        from kiube.print3d import export_print
        M = export_print(P, E, build, scale=a.escala, printer=a.impresora, outdir=str(out), paleta=a.paleta)
        report["cuerpos"] = {c["nombre"]: c["imprimible"] for c in M["cuerpos"]}
        report["estancas"] = sum(c["estanca"] for c in M["cuerpos"])
        report["n_cuerpos"] = len(M["cuerpos"])
    elif a.lod == "fem_rf":                         # modo acelerador del resonador (NGSolve): f, Q, potencia a 40 kV
        # proceso aparte: el solver necesita casi toda la memoria (7 GB) y no debe compartirla con este proceso
        import subprocess
        # la factorización queda al límite de la memoria (pico ≈ 6.4 GB): reintenta y, si no, malla algo más gruesa
        for mh in ("0.035", "0.035", "0.038"):
            r = subprocess.run([sys.executable, str(ROOT / "fem/rf_kiube.py"), "--maxh", mh, "--tag", "nominal",
                                "--figura", str(out / "modo_acelerador.png")], cwd=ROOT, capture_output=True, text=True)
            if r.returncode != -9:
                break
        print("\n".join(l for l in r.stdout.splitlines() if l.startswith("[")))
        if r.returncode != 0:
            raise SystemExit(f"fem_rf falló (código {r.returncode}): {r.stderr[-1500:]}")
        report["fem_rf"] = json.load(open(out / "resultado_nominal.json"))
        report["f_RF_objetivo_MHz"] = P["beam"]["rf_MHz"]
    elif a.lod == "stl_hd":
        from kiube.stl_hd import export_hd
        esc = f"1-{round(1 / a.escala)}" if a.escala < 1 else "1-1"
        sub = out / (f"escala_{esc}_{a.calidad}" + (f"_remallado{a.remallar:g}" if a.remallar else "")
                     + (f"_{a.paleta}" if a.paleta != "foto" else ""))
        M = export_hd(P, E, build, escala=a.escala, calidad=a.calidad, remesh_mm=a.remallar, printer=a.impresora,
                      outdir=str(sub), completo=not a.sin_completo, puente_mm=a.puente, paleta=a.paleta)
        report.update(salida=str(sub.relative_to(out)), resumen=M["resumen"], completo=M.get("completo"))
    elif a.lod == "mc":
        import shutil
        from kiube.mc import mc_volumes, export_h5m, conformal_surface_mesh, write_gdml, target_sources
        vols, info = mc_volumes(P, E, "openmc")
        report["n_volumes"] = len(vols)
        if not a.no_check:
            report["overlaps"] = overlaps(vols)
        report["invalid"] = [v["name"] for v in vols if not v["shape"].isValid()]
        # --- DAGMC (OpenMC y DAG-MCNP) + GDML desde UNA misma malla conforme (Gmsh fragment)
        tags_openmc = ["void" if v["mat"] == "vacuum" else v["mat"] for v in vols]
        vm, _ = mc_volumes(P, E, "mcnp")
        ids, tags_mcnp = {}, []
        for v in vm:
            mname = v["mat"]
            if mname == "vacuum":
                tags_mcnp.append("Vacuum"); continue
            if mname == "graveyard":
                tags_mcnp.append("Graveyard"); continue
            ids.setdefault(mname, len(ids) + 1)
            tags_mcnp.append(f"{ids[mname]}/rho:-{mats[mname]['density']}")
        ids.setdefault("air", len(ids) + 1)
        json.dump({"material_ids": ids, "implicit_complement": f"air = mat {ids['air']} (grupo mat:{ids['air']}/rho:-0.0012_comp)",
                   "densities_g_cm3": {k: mats[k]["density"] for k in ids}},
                  open(out / "dagmcnp_materials.json", "w"), indent=1)
        coords, meshes = conformal_surface_mesh(vols, h5m=[(out / "kiube_openmc.h5m", tags_openmc)])
        conformal_surface_mesh(vm, h5m=[(out / "kiube_dagmcnp.h5m", tags_mcnp)],
                               ic_tag=f"{ids['air']}/rho:-0.0012")
        ratios = {v["name"]: round(m["mesh_volume"] / v["shape"].Volume(), 5) for v, m in zip(vols, meshes)}
        report["tessellation_volume_ratio_minmax"] = [min(ratios.values()), max(ratios.values())]
        report["gdml"] = write_gdml(vols, coords, meshes, mats, out / "kiube_geant4.gdml")
        # --- fuentes, materiales y scripts
        S = target_sources(P, E)
        json.dump(json.loads(json.dumps(S, default=float)), open(out / "sources.json", "w"), indent=1)
        shutil.copy(ROOT / "materials.yaml", out / "materials.yaml")
        for f in ("run_openmc.py", "run_geant4.py", "README_MC.md"):
            shutil.copy(ROOT / "mc_templates" / f, out / f)
        shutil.copy(ROOT / "validate_g4.py", out / "validate_g4.py")
        json.dump({v["name"]: {"material": v["mat"], "volume_cm3": round(v["shape"].Volume() / 1000, 3)} for v in vols},
                  open(out / "volumes.json", "w"), indent=1)
    elif a.lod == "fem_mag":
        for tag, sector in (("sector90", True), ("half360", False)):
            doms = build_fem_mag(P, E, sector=sector)
            d = out / tag
            (d / "domains").mkdir(parents=True, exist_ok=True)
            for p in doms:
                cq.exporters.export(p["shape"], str(d / "domains" / f"{p['name']}.step"))
            save_assembly(doms, d / f"kiube_fem_mag_{tag}.step", f"KIUBE_fem_mag_{tag}")
            json.dump([{"name": p["name"], "material": p["mat"], "role": p["role"],
                        "volume_mm3": round(p["shape"].Volume(), 1)} for p in doms],
                      open(d / "domains.json", "w"), indent=1)
            report[tag] = {"n_domains": len(doms), "invalid": [p["name"] for p in doms if not p["shape"].isValid()]}
            if not a.no_check:
                report[tag]["overlaps"] = overlaps(doms)
        coil = P["coil"]
        area = coil["radial_width"] * coil["height"] * 1e-6
        (out / "SETUP_FEM_MAG.md").write_text(f"""# Configuración magnetostática (COMSOL AC/DC · Opera-3D · ANSYS Maxwell)

Unidades del CAD: mm. Mitad superior (z ≥ 0).

## Condiciones de borde
- **z = 0 (plano medio):** el campo es normal al plano (simetría par del dipolo) → H tangencial = 0.
  COMSOL: *Perfect Magnetic Conductor* (n × H = 0). Opera: condición *normal magnetic*. ANSYS: flux-normal.
- **sector90:** caras en φ = 45° y φ = 135° (centros de colina) → *Periodic Condition* (continuidad, rotación 90°).
- **Borde exterior del aire:** *Magnetic Insulation* (A×n = 0); el aire se extiende ~1.6 × el lado del yugo.

## Materiales
- Hierro (`role: iron`): curva B-H `bh_acero_bajoC_kiube` de `materials.yaml` (acero bajo C, J_s 2.12 T, CALC v0.18; reemplazar por la medida).
- Bobina (`coil_upper`): *Coil* homogeneizado, corriente azimutal.
  NI total = {coil['NI_FE_calibrado']} A·v (FE calibrado v0.18; {coil['NI_total']} por potencia) → por bobina {coil['NI_FE_calibrado']/2:.0f} A·v,
  J = {coil['NI_FE_calibrado']/2/area/1e6:.2f} A/mm² sobre {coil['radial_width']}×{coil['height']} mm.
- Aire: μr = 1. Cobre de dees, inox de cámara, Nb de blancos ≡ aire (no magnéticos, omitidos).

## Calibración sugerida
Ajustar NI hasta ⟨B⟩ = 1.335 T en r → 0 (B0 isócrono para 40.646 MHz, h = 2) y comparar
B(r, φ) con el mapa medido de SUPAG09 (B_colina ≈ 2.2 T, B_valle ≈ 0.45 T).

## Malla recomendada
- Gap de colina 24 mm: ≥ 6 elementos en z; elementos ≤ 4 mm en caras de polo, ranura e inserto.
- Bordes de polo (escalón 2 y cortes correctores): ≤ 2 mm (gradientes de campo de borde).
- Aire exterior: crecimiento 1.3–1.5.
""")
    else:
        parts, info = build(P, E, a.lod)
        report["n_parts"] = len(parts)
        report["invalid"] = [p["name"] for p in parts if not p["shape"].isValid()]
        report["masses_by_material_kg"], report["mass_by_family_kg"] = mass_report(parts, mats)
        report["targets"] = [{k: round(v, 1) if isinstance(v, float) else v for k, v in t.items()} for t in E["targets"]]
        report["ports"] = [{k: round(v, 1) if isinstance(v, float) else v for k, v in t.items()} for t in E["ports"]]
        if not a.no_check:
            report["overlaps"] = overlaps(parts)
        save_assembly(parts, out / f"kiube_{a.lod}.step", f"KIUBE_{a.lod}")
        if a.lod == "render":                      # STL por pieza + scene.json para blender/render_kiube.py
            import shutil
            from kiube.render_export import export_render
            shutil.rmtree(out / "parts", ignore_errors=True)
            report["render_piezas"] = len(export_render(P, E, build, out, variante=a.paleta)["piezas"])
        elif a.lod == "print":
            (out / "stl").mkdir(exist_ok=True)
            for p in parts:
                cq.exporters.export(p["shape"], str(out / "stl" / f"{p['name']}.stl"), tolerance=0.05, angularTolerance=0.1)
        elif a.lod == "fem_struct":
            (out / "parts").mkdir(exist_ok=True)
            for p in parts:
                cq.exporters.export(p["shape"], str(out / "parts" / f"{p['name']}.step"))
            json.dump([{"name": p["name"], "material": p["mat"], "mass_kg": round(p["mass"], 3),
                        **mats.get("fem_properties", {}).get(p["mat"], {})} for p in parts],
                      open(out / "parts_materials.json", "w"), indent=1)
        else:  # neutronics
            bymat = {}
            for p in parts:
                bymat.setdefault(p["mat"], []).append(p)
            (out / "by_material").mkdir(exist_ok=True)
            for m, ps in bymat.items():
                cq.exporters.export(cq.Compound.makeCompound([p["shape"] for p in ps]), str(out / "by_material" / f"{m}.step"))
            json.dump([{"name": p["name"], "material": p["mat"]} for p in parts], open(out / "volume_materials.json", "w"), indent=1)
            (out / "cubit_dagmc.py").write_text(f'''# Ejecutar dentro de Coreform Cubit (Tools > Play Python) desde esta carpeta.
import cubit
MATS = {json.dumps(list(bymat))}
cubit.cmd("reset"); cubit.cmd("undo off")
for m in MATS:
    before = set(cubit.get_entities("volume"))
    cubit.cmd(f'import step "by_material/{{m}}.step" heal')
    new = sorted(set(cubit.get_entities("volume")) - before)
    if new:
        cubit.cmd(f'group "mat:{{m}}" add volume ' + " ".join(map(str, new)))
cubit.cmd("imprint volume all"); cubit.cmd("merge volume all")
# graveyard (piso en z = {info["floor_z"]:.0f} mm)
cubit.cmd("brick x 6000 y 6000 z 6000"); outer = cubit.get_last_id("volume")
cubit.cmd("brick x 5800 y 5800 z 5800"); inner = cubit.get_last_id("volume")
cubit.cmd(f"subtract volume {{inner}} from volume {{outer}}")
cubit.cmd(f'group "mat:Graveyard" add volume {{cubit.get_last_id("volume")}}')
cubit.cmd('export cubit "kiube_neutronics.cub5" overwrite')
cubit.cmd('export dagmc "kiube_neutronics.h5m" faceting_tolerance 0.05 make_watertight')
''')
        json.dump({"orbit": E["orbit"][::20].round(2).tolist(), "paths": [p[::5].round(2).tolist() for p in E["paths"]],
                   "strippers": E["strippers"], "targets": E["targets"], "ports": E["ports"]},
                  open(out / "optics.json", "w"))
    json.dump(report, open(out / "report.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in report.items() if k not in ("targets", "ports")}, indent=1))


if __name__ == "__main__":
    main()
