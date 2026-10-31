"""Modelos Monte Carlo del KIUBE: DAGMC (.h5m) para OpenMC / DAG-MCNP y GDML teselado para Geant4.

Cadena: CadQuery → (vacío explícito + búnker opcional) → malla de superficie CONFORME
(imprint/fragment: las caras compartidas tienen triángulos idénticos) → h5m / GDML.
Metodología validada en Du (MIT 2024): CadQuery → cad-to-dagmc (Gmsh) → DAGMC en OpenMC,
aire como implicit complement; DAGMC ~2–3× más lento que CSG con MOAB, ≈ CSG con Embree.
"""
import math, json, pathlib, tempfile
import cadquery as cq
from .geometry import build, zcyl, mirror_z

# masas atómicas estándar (g/mol) de los elementos usados
AMASS = {"H": 1.008, "He": 4.0026, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "Na": 22.990,
         "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Ar": 39.948, "K": 39.098,
         "Ca": 40.078, "Ti": 47.867, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693,
         "Cu": 63.546, "Zn": 65.38, "Nb": 92.906, "Mo": 95.95, "Ag": 107.868, "Eu": 151.964, "W": 183.84}
ZNUM = {"H": 1, "He": 2, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Ar": 18, "K": 19, "Ca": 20, "Ti": 22, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
        "Cu": 29, "Zn": 30, "Nb": 41, "Mo": 42, "Ag": 47, "Eu": 63, "W": 74}


def mc_volumes(P, E, variant="openmc"):
    """Lista de volúmenes MC: piezas neutrónicas + vacío explícito + búnker (+ graveyard para MCNP)."""
    parts, info = build(P, E, "neutronics")
    Y, VC = P["yoke"], P["vacuum_chamber"]
    H = Y["height"] / 2
    zA = info["zA"]
    # ---- vacío: interior de la cámara + agujeros de bombeo a través de las tapas A
    vac = zcyl(VC["inner_radius"], -zA, zA)
    per = 360.0 / P["pole"]["n_sectors"]
    for k in range(P["pole"]["n_sectors"]):
        a = math.radians(P["pole"]["hill_center0"] + per / 2 + k * per)
        x, y = Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a)
        hole = zcyl(Y["pump_hole_radius"], zA, H, x, y)
        vac = vac.fuse(hole).fuse(mirror_z(hole))
    # interiores de carcasas RF y bombas difusoras (vacío real; vía Molflow se bombea desde ahí)
    X = P.get("exterior", {})
    if X:
        per_ = 360.0 / P["pole"]["n_sectors"]
        dee_az = [round(a) % 360 for a in P["dee"]["valleys"]]
        rh, hh = X["rf_housing_d"] / 2, X["rf_housing_h"]
        pr = X["pump_d"] / 2
        for k in range(P["pole"]["n_sectors"]):
            az = P["pole"]["hill_center0"] + per_ / 2 + k * per_
            a = math.radians(az)
            x, y = Y["pump_hole_r_pos"] * math.cos(a), Y["pump_hole_r_pos"] * math.sin(a)
            if round(az) % 360 in dee_az:
                cav = zcyl(rh - 20, H, H + hh - 20, x, y)
                vac = vac.fuse(cav).fuse(mirror_z(cav))
            else:
                vac = vac.fuse(zcyl(pr - 8, -H - X["pump_h"] + 15, -H, x, y))
    # perforaciones de los colimadores en la pared de la cámara: vacío hasta la ventana de Ti (camino del haz)
    from .geometry import axis_cyl
    for t in E["targets"]:
        bore = axis_cyl(P["target"]["collimator_od"] / 2, -200, 0, (t["x"], t["y"]), t["dir"])
        ro_ = VC["inner_radius"] + VC["wall"]
        vac = vac.fuse(bore.intersect(zcyl(ro_, -zA, zA)))      # solo dentro de la pared de la cámara
    vac = vac.clean()
    bb = vac.BoundingBox()
    for p in parts:
        b = p["shape"].BoundingBox()
        if b.xmax < bb.xmin or b.xmin > bb.xmax or b.ymax < bb.ymin or b.ymin > bb.ymax \
                or b.zmax < bb.zmin or b.zmin > bb.zmax:
            continue
        vac = vac.cut(p["shape"])
    vols = list(parts) + [dict(name="vacuum_interior", shape=vac, mat="vacuum", role="vacuum")]
    # ---- búnker
    V = P.get("vault", {})
    if V.get("enabled"):
        fz = info["floor_z"]
        ix, iy, ih = V["inner_x"], V["inner_y"], V["inner_height"]
        w = V["wall"]
        outer = cq.Solid.makeBox(ix + 2 * w, iy + 2 * w, ih + V["floor"] + V["ceiling"],
                                 pnt=cq.Vector(-ix / 2 - w, -iy / 2 - w, fz - V["floor"]))
        inner = cq.Solid.makeBox(ix, iy, ih, pnt=cq.Vector(-ix / 2, -iy / 2, fz))
        vols.append(dict(name="vault_concrete", shape=outer.cut(inner), mat="concrete_portland", role="shield"))
        world_bb = outer.BoundingBox()
    else:
        world_bb = None
    if variant == "mcnp":
        b = world_bb or cq.Compound.makeCompound([v["shape"] for v in vols]).BoundingBox()
        m = 100.0
        o = cq.Solid.makeBox(b.xlen + 2 * m + 400, b.ylen + 2 * m + 400, b.zlen + 2 * m + 400,
                             pnt=cq.Vector(b.xmin - m - 200, b.ymin - m - 200, b.zmin - m - 200))
        i = cq.Solid.makeBox(b.xlen + 2 * m, b.ylen + 2 * m, b.zlen + 2 * m,
                             pnt=cq.Vector(b.xmin - m, b.ymin - m, b.zmin - m))
        vols.append(dict(name="graveyard", shape=o.cut(i), mat="graveyard", role="boundary"))
    return vols, info


# ------------------------------------------------------------------- DAGMC
def export_h5m(vols, path, tags=None, min_mesh_size=0.0, max_mesh_size=25.0):
    """Escribe DAGMC en cm (scale 0.1). tags: lista de material tags (por defecto el nombre de material)."""
    from cad_to_dagmc import CadToDagmc
    m = CadToDagmc()
    tags = tags or [v["mat"] for v in vols]
    m.add_cadquery_object(cq.Compound.makeCompound([v["shape"] for v in vols]), material_tags=tags)
    return m.export_dagmc_h5m_file(filename=str(path), implicit_complement_material_tag="air",
                                   scale_factor=0.1, meshing_backend="gmsh",
                                   min_mesh_size=min_mesh_size, max_mesh_size=max_mesh_size)   # cm


# ------------------------------------------------------------------- malla conforme (gmsh)
def conformal_surface_mesh(vols, max_size=250.0, curv_n=32, threads=0, h5m=None, tight=((535.0, 585.0, 35.0), (745.0, 790.0, 40.0)), ic_tag="air"):
    """curv_n = segmentos por vuelta (flecha ≈ R·π²/(2·N²)). 'tight' = anillos (r0, r1, h_max) donde
    hay huecos estrechos: ahí el tamaño se limita para que la flecha (0.27 mm con h = 35 mm a
    R = 570) quede muy por debajo del hueco de 2 mm (sector C–cámara, cámara–bobina).
    h5m: lista de (ruta, tags) — escribe DAGMC desde ESTA misma malla (mm → cm)."""
    """Fragmenta todos los volúmenes juntos y malla superficies; devuelve nodos y triángulos
    orientados (normal hacia afuera) por volumen, en mm."""
    import gmsh
    import numpy as np
    tmp = pathlib.Path(tempfile.mkdtemp()) / "mc.step"
    cq.exporters.export(cq.Compound.makeCompound([v["shape"] for v in vols]), str(tmp))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", threads)
    dt = gmsh.model.occ.importShapes(str(tmp))
    gmsh.model.occ.synchronize()
    vol_in = [d for d in dt if d[0] == 3]
    # mapeo robusto entrada→volumen: masa + centro de masa
    ref = []
    for v in vols:
        c = v["shape"].Center() if hasattr(v["shape"], "Center") else None
        ref.append((v["shape"].Volume(), (c.x, c.y, c.z)))
    g_props = []
    for d in vol_in:
        mass = gmsh.model.occ.getMass(*d)
        cm = gmsh.model.occ.getCenterOfMass(*d)
        g_props.append((mass, cm))
    assign = {}
    used = set()
    for i, (vm, vc) in enumerate(ref):
        best, bj = None, None
        for j, (gm, gc) in enumerate(g_props):
            if j in used:
                continue
            e = abs(gm - vm) / max(vm, 1e-9) + math.dist(vc, gc) / 1000.0
            if best is None or e < best:
                best, bj = e, j
        assign[i] = vol_in[bj]
        used.add(bj)
    out, omap = gmsh.model.occ.fragment(vol_in, [])
    gmsh.model.occ.synchronize()
    in_index = {d: k for k, d in enumerate(vol_in)}
    vol_tags = {}
    for i, d in assign.items():
        vol_tags[i] = [t for (dim, t) in omap[in_index[d]] if dim == 3]
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_size)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curv_n)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    def _size(dim, tag, x, y, z, lc):
        r = math.hypot(x, y)
        for r0, r1, h in tight:
            if r0 <= r <= r1:
                return min(lc, h)
        return lc
    gmsh.model.mesh.setSizeCallback(_size)
    gmsh.model.mesh.generate(2)
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    coords = np.array(ncoord).reshape(-1, 3)
    idx = {t: k for k, t in enumerate(ntags)}
    surf_tris = {}
    for (dim, s) in gmsh.model.getEntities(2):
        et, _, en = gmsh.model.mesh.getElements(2, s)
        tri = []
        for t, nodes in zip(et, en):
            if t == 2:
                tri.extend(np.array(nodes, dtype=np.int64).reshape(-1, 3).tolist())
        surf_tris[s] = tri
    def signed_vol(tri):
        p0, p1, p2 = coords[tri[:, 0]], coords[tri[:, 1]], coords[tri[:, 2]]
        return np.einsum("ij,ij->i", p0, np.cross(p1, p2)).sum() / 6.0

    def orient_volume(surfs):
        """Orientación topológica: superficies vecinas recorren su arista común en sentidos opuestos."""
        S = [np.array([[idx[a] for a in t] for t in surf_tris[s]], dtype=np.int64).reshape(-1, 3) for s in surfs]
        edge_owner = {}
        for si, T in enumerate(S):
            cnt = {}
            for t in T:
                for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                    key = (min(a, b), max(a, b))
                    cnt[key] = cnt.get(key, 0) + 1
                    edge_owner.setdefault(key, []).append((si, a < b))
        adj = {i: [] for i in range(len(S))}
        for key, own in edge_owner.items():
            surfs_here = {}
            for si, d in own:
                surfs_here.setdefault(si, []).append(d)
            if len(surfs_here) == 2:
                (s1, d1), (s2, d2) = [(k, v[0]) for k, v in surfs_here.items()]
                same = d1 == d2
                adj[s1].append((s2, same)); adj[s2].append((s1, same))
        o = [0] * len(S)
        comps = []
        for start in range(len(S)):
            if o[start]:
                continue
            o[start] = 1; stack = [start]; comp = [start]
            while stack:
                u = stack.pop()
                for v_, same in adj[u]:
                    want = -o[u] if same else o[u]
                    if not o[v_]:
                        o[v_] = want; stack.append(v_); comp.append(v_)
            comps.append(comp)
        oriented = [T if o[i] > 0 else T[:, [0, 2, 1]] for i, T in enumerate(S)]
        # signo por componente: la de mayor |V| positiva (envolvente), las demás negativas (cavidades)
        cv = [signed_vol(np.vstack([oriented[i] for i in c])) if c else 0.0 for c in comps]
        big = int(np.argmax(np.abs(cv))) if cv else 0
        for ci, c in enumerate(comps):
            want_pos = ci == big
            if (cv[ci] > 0) != want_pos:
                for i in c:
                    oriented[i] = oriented[i][:, [0, 2, 1]]
        return np.vstack(oriented) if oriented else np.zeros((0, 3), dtype=np.int64), len(comps)

    def surf_flip(s, sign):
        """True si los triángulos de la superficie s deben invertirse para apuntar hacia afuera.
        Normal geométrica de la cara CAD (gmsh.getNormal) × signo de frontera del volumen."""
        T = surf_tris[s]
        k = len(T) // 2
        p = np.array([coords[idx[a]] for a in T[k]])
        nm = np.cross(p[1] - p[0], p[2] - p[0])
        c = p.mean(0)
        uv = gmsh.model.getParametrization(2, s, c.tolist())
        ng = np.array(gmsh.model.getNormal(s, uv))
        return float(np.dot(nm, ng)) * sign < 0

    result = []
    for i in range(len(vols)):
        tri_l = []
        for vt in vol_tags[i]:
            for (d, s_) in gmsh.model.getBoundary([(3, vt)], combined=False, oriented=True, recursive=False):
                st = surf_tris[abs(s_)]
                if not st:
                    continue
                f = surf_flip(abs(s_), 1 if s_ > 0 else -1)
                tri_l.extend([[t[0], t[2], t[1]] for t in st] if f else st)
        tri = np.array([[idx[a] for a in t] for t in tri_l], dtype=np.int64).reshape(-1, 3)
        sv = signed_vol(tri) if len(tri) else 0.0
        if sv < 0:                         # convención de signo global opuesta → invertir todo
            tri = tri[:, [0, 2, 1]]; sv = -sv
        # verificación: aristas dirigidas repetidas = orientación inconsistente
        from collections import Counter
        ed = Counter()
        for t in tri:
            ed[(t[0], t[1])] += 1; ed[(t[1], t[2])] += 1; ed[(t[2], t[0])] += 1
        bad = sum(1 for c in ed.values() if c > 1)
        if bad:                            # respaldo topológico
            surfs = []
            for vt in vol_tags[i]:
                surfs += [abs(s_) for (d, s_) in gmsh.model.getBoundary([(3, vt)], combined=False, oriented=False,
                                                                        recursive=False) if surf_tris[abs(s_)]]
            tri2, _ = orient_volume(surfs)
            ed2 = Counter()
            for t in tri2:
                ed2[(t[0], t[1])] += 1; ed2[(t[1], t[2])] += 1; ed2[(t[2], t[0])] += 1
            bad2 = sum(1 for c in ed2.values() if c > 1)
            if bad2 < bad:
                tri, bad, sv = tri2, bad2, signed_vol(tri2)
        result.append(dict(tri=tri, mesh_volume=sv, bad_edges=bad))

    if h5m:
        from cad_to_dagmc.core import export_gmsh_object_to_dagmc_h5m_file
        import h5py
        order = [t for (d, t) in gmsh.model.getEntities(3)]
        owner = {}
        for i, tl in vol_tags.items():
            for t in tl:
                owner[t] = i
        for path, tags in h5m:
            mt = [tags[owner[t]] for t in order]
            export_gmsh_object_to_dagmc_h5m_file(material_tags=mt, implicit_complement_material_tag=ic_tag,
                                                 filename=str(path))
            with h5py.File(str(path), "r+") as f:
                f["tstt/nodes/coordinates"][...] = f["tstt/nodes/coordinates"][...] * 0.1   # mm → cm
    gmsh.finalize()
    return coords, result


# ------------------------------------------------------------------- GDML (Geant4)
def _material_xml(name, m):
    """Devuelve (elementos/isótopos xml, material xml) para un material de materials.yaml."""
    defs, mat = [], []
    if name == "vacuum":
        return [], ['<material name="vacuum" Z="1"><D value="1e-25" unit="g/cm3"/><atom value="1.008"/></material>']
    rho = m["density"]
    if "isotopic" in m:              # agua: H2 + O (natural o enriquecido)
        iso = m["isotopic"]
        o16, o18 = iso.get("O16", 0.0), iso.get("O18", 0.0)
        en = f"O_{name}"
        defs += ['<isotope name="O16" Z="8" N="16"><atom value="15.9949"/></isotope>',
                 '<isotope name="O18" Z="8" N="18"><atom value="17.9992"/></isotope>',
                 f'<element name="{en}" formula="O"><fraction n="{o16 / (o16 + o18):.6f}" ref="O16"/>'
                 f'<fraction n="{o18 / (o16 + o18):.6f}" ref="O18"/></element>']
        mat.append(f'<material name="{name}" state="liquid"><D value="{rho}" unit="g/cm3"/>'
                   f'<composite n="2" ref="el_H"/><composite n="1" ref="{en}"/></material>')
        return defs, mat
    wt = m["wt"]
    tot = sum(wt.values())
    state = "gas" if rho < 0.05 else "solid"
    fr = "".join(f'<fraction n="{v / tot:.8f}" ref="el_{k}"/>' for k, v in wt.items() if v > 0)
    mat.append(f'<material name="{name}" state="{state}"><D value="{rho}" unit="g/cm3"/>{fr}</material>')
    return defs, mat


def write_gdml(vols, coords, meshes, mats, path, world_margin=500.0):
    import numpy as np
    used = sorted({v["mat"] for v in vols if v["mat"] not in ("graveyard",)} | {"air"})
    els = sorted({e for k in used for e in mats.get(k, {}).get("wt", {})} | {"H", "N", "O", "Ar"})
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
             'xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">',
             "<!-- Cyclone KIUBE 18 MeV — modelo paramétrico público (v0.3). Unidades mm. "
             "Teselado conforme (Gmsh fragment): caras compartidas con triángulos idénticos. -->"]
    # define: nodos usados
    usedn = np.unique(np.concatenate([m["tri"].ravel() for m in meshes if len(m["tri"])]))
    lines.append("<define>")
    for n in usedn:
        x, y, z = coords[n]
        lines.append(f'<position name="n{n}" unit="mm" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}"/>')
    lines.append("</define>")
    # materiales
    lines.append("<materials>")
    for e in els:
        lines.append(f'<element name="el_{e}" formula="{e}" Z="{ZNUM[e]}"><atom value="{AMASS[e]}"/></element>')
    for k in used:
        d, m = _material_xml(k, mats[k])
        lines += d + m
    lines.append("</materials>")
    # sólidos
    lines.append("<solids>")
    allp = coords[usedn]
    half = np.abs(allp).max(0) + world_margin          # mundo centrado en el origen del CAD
    lines.append(f'<box name="world_s" lunit="mm" x="{2*half[0]:.3f}" y="{2*half[1]:.3f}" z="{2*half[2]:.3f}"/>')
    for v, m in zip(vols, meshes):
        if v["mat"] == "graveyard" or not len(m["tri"]):
            continue
        lines.append(f'<tessellated name="{v["name"]}_s" aunit="deg" lunit="mm">')
        for a, b, cc in m["tri"]:
            lines.append(f'<triangular vertex1="n{a}" vertex2="n{b}" vertex3="n{cc}" type="ABSOLUTE"/>')
        lines.append("</tessellated>")
    lines.append("</solids>")
    # estructura
    lines.append("<structure>")
    for v in vols:
        if v["mat"] == "graveyard":
            continue
        lines.append(f'<volume name="{v["name"]}"><materialref ref="{v["mat"]}"/><solidref ref="{v["name"]}_s"/></volume>')
    lines.append('<volume name="World"><materialref ref="air"/><solidref ref="world_s"/>')
    for v in vols:
        if v["mat"] == "graveyard":
            continue
        lines.append(f'<physvol name="{v["name"]}_pv"><volumeref ref="{v["name"]}"/></physvol>')
    lines.append("</volume></structure>")
    lines.append('<setup name="Default" version="1.0"><world ref="World"/></setup></gdml>')
    pathlib.Path(path).write_text("\n".join(lines))
    return dict(world_half_mm=half.tolist(), n_nodes=int(len(usedn)),
                n_triangles=int(sum(len(m["tri"]) for m in meshes)))


# ------------------------------------------------------------------- término fuente
def target_sources(P, E, current_uA=None, dual=None, normalizacion="BON21"):
    """Fuentes de neutrones en el centro de la cavidad de cada blanco.
    Espectro: H2-18O, 18 MeV p, Bonvin et al. 2021 Fig.7 (FLUKA), extraído del PDF vectorial.
    Rendimiento integrado 4.9e-3 n/p. Bonvin halló que FLUKA subestima la activación ~2.5×
    (emisión hacia atrás subestimada, Bakhtiari 2020): factor disponible como 'bonvin_correction'."""
    from .geometry import axis_cyl  # noqa: F401  (misma convención de ejes)
    spec = json.load(open(pathlib.Path(__file__).resolve().parent.parent / "data" / "h2o18_spectrum_bon21.json"))
    from .geometry import target_cavity
    TG = P["target"]
    PR = P.get("proyecto", {})
    current_uA = current_uA if current_uA is not None else PR.get("I_max_uA", 150.0)    # PROP: KIUBE 180
    dual = dual if dual is not None else PR.get("doble_haz", True)
    # normalización del rendimiento n/p: 'BON21' (FLUKA, espectro) o 'fuente_rad' (params.fuente_rad.n_por_p_18O)
    ypp = spec["yield_n_per_proton"] if normalizacion == "BON21" else P["fuente_rad"]["n_por_p_18O"]   # 'fuente_rad'
    t1 = TG["vac_window_um"] / 1000
    s3 = t1 + TG["he_gap"] + TG["tgt_window_um"] / 1000
    pps = current_uA * 1e-6 / 1.602176634e-19
    n_active = 2 if dual else 1                     # KIUBE: 2 blancos simultáneos (doble extracción)
    out = []
    for i, t in enumerate(E["targets"]):
        ttype = TG["types"][i] if i < len(TG["types"]) else "F18"
        a = math.radians(t["dir"])
        zc = s3 + target_cavity(TG, ttype)[1] / 2
        pos = [t["x"] + zc * math.cos(a), t["y"] + zc * math.sin(a), 0.0]
        out.append(dict(name=f"T{i+1}", type=ttype, pos_mm=pos, beam_dir_deg=t["dir"],
                        protons_per_s_if_active=pps / n_active,
                        neutrons_per_s_if_active=pps / n_active * ypp
                        if ttype == "F18" else None))
    return dict(sources=out, spectrum=spec, current_uA=current_uA, simultaneous_targets=n_active,
                normalizacion=normalizacion, yield_n_per_proton=ypp,
                bonvin_correction=2.5, note="activar un subconjunto de blancos según el escenario")
