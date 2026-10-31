"""Vacío del Cyclone KIUBE: exportación a Molflow+ y pérdidas de H⁻ por stripping en gas residual.

Cadena:  CAD (vacío explícito) → superficies clasificadas por su vecino (bombas, acero, cobre,
         inoxidable, fuentes de H₂) → Molflow+ (presión 3D)  →  pérdidas de H⁻ a lo largo de la
         espiral de aceleración  →  puntos de impacto de H⁰ en la pared  →  término fuente "BL"
         (pérdidas de haz) para Geant4.

Mientras no haya resultado de Molflow se usa un estimador analítico de presión media (Q/S_eff)
calibrado con la presión base medida (TUD03: 3·10⁻⁷ mbar con RF) y una sección eficaz con la
forma σ ∝ 1/β² (Born) cuya MAGNITUD se calibra con la fracción de pérdida por gas declarada
(por defecto 15 %: transmisión KIUBE 81 % [TUD03], AMIT 14.6 % [arXiv:2103.15583]).
Nada de esto sustituye una tabla σ(E) medida (Nakai/Tabata-Shirai) ni el campo de presión de Molflow.
"""
import json, math, pathlib, re
import numpy as np

KB = 1.380649e-23
R_GAS = 8.314462618


def vbar_cm_s(M_g_mol=2.016, T=293.15):
    return math.sqrt(8 * R_GAS * T / (math.pi * M_g_mol * 1e-3)) * 100.0


# ------------------------------------------------------------------ superficies del vacío
def classify_vacuum_surfaces(vols, P, curv_n=24, max_size=60.0):
    """Malla conforme y clasificación de cada superficie del volumen de vacío por su vecino.
    Devuelve grupos {nombre: triángulos (N×3×3) en mm con normal HACIA el gas} + metadatos."""
    import gmsh, cadquery as cq, tempfile
    global P_H, P_PUMP, P_RFH
    P_H, P_PUMP, P_RFH = P["yoke"]["height"] / 2, P["exterior"]["pump_h"], P["exterior"]["rf_housing_h"]
    iv = next(i for i, v in enumerate(vols) if v["mat"] == "vacuum")
    vb = vols[iv]["shape"].BoundingBox()
    near = [i for i, v in enumerate(vols) if i == iv or not (
        v["shape"].BoundingBox().xmax < vb.xmin or v["shape"].BoundingBox().xmin > vb.xmax or
        v["shape"].BoundingBox().ymax < vb.ymin or v["shape"].BoundingBox().ymin > vb.ymax or
        v["shape"].BoundingBox().zmax < vb.zmin or v["shape"].BoundingBox().zmin > vb.zmax)]
    sub = [vols[i] for i in near]
    tmp = pathlib.Path(tempfile.mkdtemp()) / "vac.step"
    cq.exporters.export(cq.Compound.makeCompound([v["shape"] for v in sub]), str(tmp))
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    dt = [d for d in gmsh.model.occ.importShapes(str(tmp)) if d[0] == 3]
    gmsh.model.occ.synchronize()
    out, omap = gmsh.model.occ.fragment(dt, [])
    gmsh.model.occ.synchronize()
    owner = {}
    for k, m in enumerate(omap):
        for d, t in m:
            if d == 3:
                owner[t] = sub[k]
    vac_tags = [t for t, v in owner.items() if v["mat"] == "vacuum"]
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_size)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curv_n)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    nt, nc, _ = gmsh.model.mesh.getNodes()
    C = np.array(nc).reshape(-1, 3); idx = {t: k for k, t in enumerate(nt)}
    groups, meta = {}, {}
    for vt in vac_tags:
        for (d, s) in gmsh.model.getBoundary([(3, vt)], combined=False, oriented=False, recursive=False):
            s = abs(s)
            up, _ = gmsh.model.getAdjacencies(2, s)
            if len(up) == 2 and all(u in vac_tags for u in up):
                continue                                   # cara interna entre fragmentos del vacío
            others = [owner[u] for u in up if u not in vac_tags]
            et, _, en = gmsh.model.mesh.getElements(2, s)
            tri = np.array(en[0], dtype=np.int64).reshape(-1, 3) if len(en) else np.zeros((0, 3), int)
            if not len(tri):
                continue
            T = C[np.vectorize(idx.get)(tri)]
            cen = T.reshape(-1, 3).mean(0)
            if others:
                o = others[0]
                g = {"steel_1010": "pared_acero", "copper_ofe": "pared_cobre", "copper_dee_caja": "pared_cobre", "stainless_304L": "pared_inox",
                     "aluminum_6061": "pared_aluminio", "titanium": "ventana_ti", "graphite_foil": "lamina_c",
                     "niobium": "pared_niobio"}.get(o["mat"], "pared_otro")
                az = int(round((math.degrees(math.atan2(cen[1], cen[0])) % 360) / 90.0) * 90) % 360
                m_is = re.match(r"^ion_source_(\d)$", o["name"])
                if m_is and math.hypot(cen[0], cen[1]) < 60:
                    g = "fuente_h2_" + m_is.group(1)          # chimenea de la fuente PIG (rendija de H₂)
                elif o["name"].startswith("diffusion_pump") and cen[2] < -P_H - P_PUMP + 30:
                    g = f"bomba_dp_az{az:03d}"                      # entrada de la difusora
                elif o["name"].startswith("rf_housing") and cen[2] > P_H + P_RFH - 40:
                    g = f"bomba_rf_az{az:03d}"                      # bomba montada en la carcasa RF (EST)
            else:
                r = math.hypot(cen[0], cen[1]); az = math.degrees(math.atan2(cen[1], cen[0])) % 360
                if r < 60:
                    g = "holgura_fuentes"               # anillo alrededor de las fuentes (termina en el posicionador)
                else:
                    g = f"boca_z{'+' if cen[2] > 0 else '-'}_az{int(round(az / 90.0) * 90) % 360:03d}"
            # orientar hacia el gas: el centroide desplazado ε por la normal debe quedar dentro del vacío
            n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
            k0 = int(np.argmax(np.linalg.norm(n, axis=1)))
            c0 = T[k0].mean(0); n0 = n[k0] / np.linalg.norm(n[k0])
            inside = any(gmsh.model.isInside(3, vt, (c0 + 0.05 * n0).tolist()) for vt in vac_tags)
            if not inside:
                T = T[:, [0, 2, 1], :]
            groups.setdefault(g, []).append(T)
            meta.setdefault(g, {"area_cm2": 0.0})["area_cm2"] += float(np.linalg.norm(n, axis=1).sum() / 2 / 100)
    gmsh.finalize()
    return {g: np.vstack(v) for g, v in groups.items()}, meta


def write_stl_ascii(path, tris_cm, name):
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for t in tris_cm:
            n = np.cross(t[1] - t[0], t[2] - t[0]); ln = np.linalg.norm(n); n = n / ln if ln > 0 else n
            f.write(f"facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n outer loop\n")
            for p in t:
                f.write(f"  vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}\n")
            f.write(" endloop\nendfacet\n")
        f.write(f"endsolid {name}\n")


def export_molflow(P, groups, meta, outdir):
    """STL único (cm, triángulos ordenados por grupo, normales reorientadas HACIA EL GAS por cáscara
    cerrada) + mapa de índices + un archivo de parámetros MolflowCLI por especie:
      kiube_molflow_H2O_10h.txt — desgasificación de metales a t = 10 h (CERN-CH Ec.13), bombas para H₂O
      kiube_molflow_H2_fuente.txt — gas de la fuente PIG activa, bombas para H₂."""
    import shutil
    from .tpmc import load_gas_mesh, sticking_for_speed
    V = P["vacuum"]
    out = pathlib.Path(outdir)
    shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True, exist_ok=True)
    order = sorted(groups)
    allt, ranges, i0 = [], {}, 1
    for g in order:
        t = groups[g] * 0.1                                   # mm → cm (unidad interna de Molflow)
        allt.append(t)
        ranges[g] = [i0, i0 + len(t) - 1]
        i0 += len(t)
    write_stl_ascii(out / "kiube_vacio_cm.stl", np.vstack(allt), "kiube_vacuum")
    m = load_gas_mesh(out / "kiube_vacio_cm.stl")             # normales al gas sin reordenar facetas
    T_or = m.vertices[m.faces]
    write_stl_ascii(out / "kiube_vacio_cm.stl", T_or, "kiube_vacuum")
    for g, (a, b) in ranges.items():
        write_stl_ascii(out / f"grupo_{g}.stl", T_or[a - 1:b], g)
    vb = vbar_cm_s(V["gas_mass_g_mol"], V["temperature_K"])
    pumps = {g: ranges[g] for g in ranges if g in V["active_pumps"]}
    q10 = V["outgassing_q1h_H2O"] * 10.0 ** (-V["outgassing_alpha"])
    head = ["# MolflowCLI --setParamsByFile — KIUBE (generado). Índices de facetas 1-based en el orden del STL.",
            "# Importar kiube_vacio_cm.stl en Molflow (File > Import > STL, unidad cm) SIN 'Collapse'; guardar .xml;",
            "# luego: molflowCLI -f kiube.xml --setParamsByFile <este archivo> -t 600 --writeFacetQuantities"]
    LH2O = head + ["# especie: H2O; desgasificación de metales sin hornear a t = 10 h [Chiggiato, arXiv:2006.07124, Ec.13]",
                   "simulation.mass=18.015"]
    LH2 = head + ["# especie: H2; gas de la fuente PIG activa", "simulation.mass=2.016"]
    info = {"grupos": {}, "v_media_cm_s": vb}
    for g, (a, b) in ranges.items():
        A = meta[g]["area_cm2"]
        e = {"facetas": [a, b], "area_cm2": round(A, 1)}
        if g in pumps:
            s2 = sticking_for_speed(V["pump_speed_l_s"], A, 2.016, V["temperature_K"])
            s18 = sticking_for_speed(V["pump_speed_H2O_l_s"], A, 18.015, V["temperature_K"])
            LH2.append(f"facet.{a}-{b}.sticking={s2:.5f}"); LH2O.append(f"facet.{a}-{b}.sticking={s18:.5f}")
            e.update(bomba=True, velocidad_l_s={"H2": V["pump_speed_l_s"], "H2O": V["pump_speed_H2O_l_s"]},
                     sticking={"H2": round(s2, 5), "H2O": round(s18, 5)})
        elif g.startswith(("boca_", "bomba_", "holgura_")):
            for L in (LH2, LH2O):
                L.append(f"facet.{a}-{b}.sticking=0")
            e.update(bomba=False, nota="boca cerrada (brida ciega / carcasa RF sin bomba)")
        elif g.startswith("fuente_h2"):
            Q = V["h2_sccm_per_active_source"] * 0.01689 if g.endswith(V["active_source"]) else 0.0
            LH2.append(f"facet.{a}-{b}.outgassing={Q / max(1, b - a + 1):.6e}")
            e.update(carga_h2_mbar_l_s=Q)
        else:
            q = q10 * V["outgassing_factor"].get(g, 1.0)
            LH2O.append(f"facet.{a}-{b}.specificOutgassing={q:.3e}")
            e.update(desgasificacion_H2O_10h_mbar_l_s_cm2=q)
        info["grupos"][g] = e
    (out / "kiube_molflow_H2O_10h.txt").write_text("\n".join(LH2O) + "\n")
    (out / "kiube_molflow_H2_fuente.txt").write_text("\n".join(LH2) + "\n")
    json.dump(info, open(out / "kiube_molflow_grupos.json", "w"), indent=1, ensure_ascii=False)
    # geometría lista para Molflow 2.11 (sin pasar por File > Import > STL): mismo orden de facetas que el STL,
    # así que los archivos --setParamsByFile valen tal cual. Se guarda comprimida (molflowCLI lee .zip).
    import zipfile
    from .molflow import write_xml
    xml = out / "kiube_vacio.xml"
    write_xml(xml, T_or, mass=V["gas_mass_g_mol"], temperature=V["temperature_K"], structure_name="kiube_vacio")
    with zipfile.ZipFile(out / "kiube_vacio.zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.write(xml, "kiube_vacio.xml")
    xml.unlink()
    info["xml_molflow"] = "kiube_vacio.zip"
    return info


def analytic_pressure(P, info):
    """Estimación rápida Q/S sin geometría (la presión de referencia es la del TPMC)."""
    V = P["vacuum"]
    vb = info["v_media_cm_s"]
    Q_out = sum(g["area_cm2"] * g.get("desgasificacion_H2O_10h_mbar_l_s_cm2", 0.0) for g in info["grupos"].values())
    Q_h2 = sum(g.get("carga_h2_mbar_l_s", 0.0) for g in info["grupos"].values())
    S2 = S18 = 0.0
    for name, g in info["grupos"].items():
        if g.get("bomba"):
            C2 = g["area_cm2"] * vb / 4 / 1000.0 * V["port_transmission"]
            C18 = C2 * math.sqrt(2.016 / 18.015)
            S2 += 1.0 / (1.0 / V["pump_speed_l_s"] + 1.0 / C2)
            S18 += 1.0 / (1.0 / V["pump_speed_H2O_l_s"] + 1.0 / C18)
    return dict(Q_H2O_10h_mbar_l_s=Q_out, Q_h2_mbar_l_s=Q_h2, S_eff_H2_l_s=S2, S_eff_H2O_l_s=S18,
                P_H2O_10h_mbar=Q_out / S18 if S18 else float("inf"),
                P_H2_operacion_mbar=Q_h2 / S2 if S2 else float("inf"))


# ------------------------------------------------------------------ pérdidas de H⁻
def lorentz_lifetime_s(T_MeV, B_T, m=939.294):
    """Vida media en el sistema propio (Keating 1995; constantes de arXiv:2103.15583, ec. 20)."""
    g = 1 + T_MeV / m; b = math.sqrt(1 - 1 / g ** 2)
    Ef = g * b * 2.99792458e8 * B_T
    A1, A2 = 2.714e-6, 4.47407e9
    return (A1 / Ef) * math.exp(A2 / Ef), Ef


def hminus_losses(P, E, pressure_mbar, loss_fraction=None, n_samples=20000, seed=7):
    """Espiral de aceleración: dN/N = −n σ(E) dl, σ = σ_ref (β_ref/β)², calibrada a loss_fraction.
    Devuelve pérdida por energía y puntos de impacto de H⁰ (trayectoria tangente recta) en la pared."""
    V, b = P["vacuum"], P["beam"]
    m = 939.294
    dE = b["dE_turn_keV"] / 1000.0
    T0 = V["injection_energy_MeV"]
    Ts = np.arange(T0, b["T_MeV"], dE)
    p18 = math.sqrt(b["T_MeV"] ** 2 + 2 * b["T_MeV"] * m)
    orb = np.asarray(E["orbit"]) ; r18 = np.hypot(orb[:, 0], orb[:, 1]); th = np.arctan2(orb[:, 1], orb[:, 0])
    L18 = float(np.sum(np.hypot(*np.diff(np.vstack([orb, orb[:1]]), axis=0).T)))   # mm
    p = np.sqrt(Ts ** 2 + 2 * Ts * m); g = 1 + Ts / m; beta = np.sqrt(1 - 1 / g ** 2)
    Lturn_cm = L18 * (p / p18) / 10.0
    r_turn = float(r18.mean()) * (p / p18)                       # radio medio de cada vuelta, mm
    Pv = np.array([pressure_mbar(r) for r in r_turn]) if callable(pressure_mbar) else np.full(len(Ts), pressure_mbar)
    n = Pv * 100.0 / (KB * 293.15) * 1e-6 * 2                     # átomos/cm³ (H₂ → 2 átomos)
    shape = (beta[-1] / beta) ** 2                                # σ/σ_ref, σ_ref en 18 MeV
    tau = n * shape * Lturn_cm                                    # por σ_ref
    f = loss_fraction if loss_fraction is not None else V["gas_loss_fraction"]
    sig_ref = -math.log(1 - f) / tau.sum()                        # cm², calibrada
    surv = np.exp(-np.cumsum(tau * sig_ref))
    dloss = np.concatenate([[1 - surv[0]], surv[:-1] - surv[1:]])
    # muestreo de puntos de pérdida y vuelo recto tangente del H⁰ hasta r = pared interior de la cámara
    rng = np.random.default_rng(seed)
    Rw = P["vacuum_chamber"]["inner_radius"]
    ks = rng.choice(len(Ts), size=n_samples, p=dloss / dloss.sum())
    js = rng.integers(0, len(orb), size=n_samples)
    hits = []
    for k, j in zip(ks, js):
        sc = p[k] / p18
        x, y = orb[j] * sc
        j2 = (j + 1) % len(orb)
        tx, ty = (orb[j2] - orb[j]); tn = math.hypot(tx, ty); tx, ty = tx / tn, ty / tn
        bq = x * tx + y * ty
        t = -bq + math.sqrt(bq * bq - (x * x + y * y - Rw * Rw))
        hx, hy = x + t * tx, y + t * ty
        cosi = abs((hx * tx + hy * ty) / Rw)
        hits.append((float(Ts[k]), hx, hy, cosi))
    hits = np.array(hits)
    return dict(sigma_ref_18MeV_cm2=sig_ref, sigma_1MeV_cm2=sig_ref * float(((beta[-1]) / math.sqrt(1 - 1 / (1 + 1 / m) ** 2)) ** 2),
                fraccion_perdida=float(1 - surv[-1]), vueltas=len(Ts),
                perdida_por_MeV=np.histogram(Ts, bins=np.arange(0, 19, 1.0), weights=dloss)[0].tolist(),
                energia_media_impacto_MeV=float(hits[:, 0].mean()),
                fraccion_impactos_E_gt_5MeV=float((hits[:, 0] > 5).mean()),
                impactos=hits)


# ------------------------------------------------------------------ resultados de Molflow → P(r)
def pressure_profile_from_molflow(csv_path, grupos_json, stl_cm, z_band_mm=40.0, bins_mm=25.0):
    """Lee el CSV de MolflowCLI (--writeFacetQuantities o --writeFacetDetails) y devuelve P(r) en mbar.
    Toma las facetas de pared cercanas al plano medio (|z| < z_band: caras de polo, dees) como
    muestra de la presión en el camino del haz. Columnas detectadas por nombre (contienen 'pressure'
    y 'facet'/'#'); el índice de faceta es 1-based y sigue el orden del STL exportado."""
    import csv, trimesh
    rows = list(csv.reader(open(csv_path, newline="")))
    hdr = [h.strip().lower() for h in rows[0]]
    ip = next(i for i, h in enumerate(hdr) if "pressure" in h)
    ifc = next((i for i, h in enumerate(hdr) if "facet" in h or h in ("#", "id", "index")), 0)
    pres = {}
    for r in rows[1:]:
        try:
            pres[int(float(r[ifc]))] = float(r[ip])
        except (ValueError, IndexError):
            continue
    m = trimesh.load(stl_cm, process=False)
    cen = m.triangles_center * 10.0                       # cm → mm
    info = json.load(open(grupos_json))["grupos"]
    walls = [(a, b) for g, v in info.items() if g.startswith("pared_") for a, b in [v["facetas"]]]
    rr, pp = [], []
    for a, b in walls:
        for f in range(a, b + 1):
            c = cen[f - 1]
            if abs(c[2]) < z_band_mm and f in pres:
                rr.append(math.hypot(c[0], c[1])); pp.append(pres[f])
    rr, pp = np.array(rr), np.array(pp)
    edges = np.arange(0, rr.max() + bins_mm, bins_mm)
    prof = [(0.5 * (edges[k] + edges[k + 1]), float(pp[(rr >= edges[k]) & (rr < edges[k + 1])].mean()))
            for k in range(len(edges) - 1) if ((rr >= edges[k]) & (rr < edges[k + 1])).any()]
    R = np.array([p[0] for p in prof]); Pm = np.array([p[1] for p in prof])
    return lambda r_mm: float(np.interp(r_mm, R, Pm)), prof
