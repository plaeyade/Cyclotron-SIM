"""Integración con Molflow+ 2.11 (CERN): escritura del XML de geometría, ejecución de molflowCLI y lectura
de resultados. Esquema obtenido de Molflow 2.11.1 (el XML que el propio molflowCLI guarda, «version 21101»).

Unidades del XML de Molflow (verificadas con casos analíticos en validation/molflow_vs_tpmc.py):
  - coordenadas en cm;
  - Outgassing constValue: Pa·m³/s POR FACETA (1 Pa·m³/s = 10 mbar·l/s); con un parámetro temporal, el valor
    del parámetro multiplica a constValue;
  - Sticking: probabilidad [0, 1]; masa del gas en g/mol; temperatura en K.
  - desType: 0 ninguna, 1 uniforme, 2 coseno (ley de Knudsen, la física correcta), 3 coseno^N.
Presión por faceta desde los resultados: P[Pa] = Σv_ort · (moléculas por partícula de prueba) · M/6e26 / A[m²]
(misma constante de masa que usa Molflow); en momentos dividir además por la ventana temporal.

Uso típico:
    from kiube.molflow import write_xml, run_cli, read_results
    write_xml("caso.xml", tris_cm, sticking=s, outgassing_pa_m3_s=q, mass=18.015)
    out = run_cli("caso.xml", ndes=1e6, outdir="res")
    R = read_results(out)              # R["moments"][0]["P_mbar"] → array por faceta (flujo constante)
"""
from __future__ import annotations

import os, pathlib, shutil, subprocess
import xml.etree.ElementTree as ET
import numpy as np

AMU_MOLFLOW = 1.0 / 6e26          # kg por (g/mol) que usa Molflow en la presión (≈ 1.6667e-27)


def find_cli():
    """molflowCLI: variable MOLFLOW_CLI, PATH, o la carpeta molflow/ del proyecto o junto a él."""
    cand = [os.environ.get("MOLFLOW_CLI"), shutil.which("molflowCLI")]
    root = pathlib.Path(__file__).resolve().parents[1]
    for d in [root / "molflow", root.parent / "molflow_pkg"]:
        if d.exists():
            cand += [str(p) for p in d.rglob("molflowCLI")]
    for c in cand:
        if c and pathlib.Path(c).exists():
            return str(c)
    return None


def _merge_vertices(tris, decimals=7):
    P = np.asarray(tris, float).reshape(-1, 3)
    key = np.round(P, decimals)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    return uniq, inv.reshape(-1, 3)


def write_xml(path, tris_cm, sticking=None, outgassing_pa_m3_s=None, mass=28.0, temperature=293.15,
              opacity=None, two_sided=None, outgassing_param=None, params=None, moments=None, time_window=1e-10,
              const_flow=True, structure_name="kiube"):
    """tris_cm: (N,3,3) triángulos en cm (normal = lado del gas, como en Molflow).
    params: {nombre: [(t_s, valor), ...]} → parámetros temporales; outgassing_param: índice de parámetro por faceta
    (−1 = constante). moments: lista de cadenas «t0,dt,t1» (UserMoments de Molflow)."""
    V, F = _merge_vertices(tris_cm)
    n = len(F)
    s = np.zeros(n) if sticking is None else np.broadcast_to(np.asarray(sticking, float), (n,))
    q = np.zeros(n) if outgassing_pa_m3_s is None else np.broadcast_to(np.asarray(outgassing_pa_m3_s, float), (n,))
    op = np.ones(n) if opacity is None else np.broadcast_to(np.asarray(opacity, float), (n,))
    ts = np.zeros(n, int) if two_sided is None else np.broadcast_to(np.asarray(two_sided, int), (n,))
    pid = -np.ones(n, int) if outgassing_param is None else np.broadcast_to(np.asarray(outgassing_param, int), (n,))
    out = ['<?xml version="1.0"?>', '<SimulationEnvironment type="MolFlow" version="21101">', "<Geometry>",
           f'<Vertices nb="{len(V)}">']
    out += [f'<Vertex id="{i}" x="{v[0]:.9g}" y="{v[1]:.9g}" z="{v[2]:.9g}" />' for i, v in enumerate(V)]
    out += ["</Vertices>", f'<Structures nb="1"><Structure id="0" name="{structure_name}" /></Structures>',
            f'<Facets nb="{n}">']
    for i in range(n):
        a, b, c = F[i]
        out.append(
            f'<Facet id="{i}"><Sticking constValue="{s[i]:.9g}" parameterId="-1" />'
            f'<Opacity constValue="{op[i]:.9g}" parameterId="-1" is2sided="{ts[i]}" />'
            f'<Outgassing constValue="{q[i]:.9g}" parameterId="{pid[i]}" desType="2" desExponent="0" '
            f'outgassingPerAreaMode="false" hasOutgassingFile="0" useOutgassingFile="0" delayOutg="false" '
            f'initialParticles="0" initialParticlesPerAreaMode="false" />'
            f'<Temperature value="{temperature}" accFactor="1" accGovernsReflection="0" />'
            f'<Reflection diffusePart="1" specularPart="0" cosineExponent="0" enableSojournTime="0" '
            f'sojournFreq="10000000000000" sojournE="100" />'
            f'<Structure inStructure="0" linksTo="0" /><Teleport target="0" /><Motion isMoving="0" />'
            f'<Recordings><Profile type="0" numSlices="100" name="none" /><Texture hasMesh="false" texDimX="0" '
            f'texDimY="0" countDes="0" countAbs="0" countRefl="0" countTrans="0" countDir="0" countAC="0" />'
            f'</Recordings><ViewSettings textureVisible="true" volumeVisible="true" />'
            f'<Indices nb="3"><Indice id="0" vertex="{a}" /><Indice id="1" vertex="{b}" /><Indice id="2" vertex="{c}" />'
            f'</Indices></Facet>')
    out += ["</Facets>", "</Geometry>", "<MolflowSimuSettings>",
            f'<Gas mass="{mass}" enableDecay="0" halfLife="1" />',
            f'<TimeSettings timeWindow="{time_window}" useMaxwellDistr="1" calcConstFlow="{int(const_flow)}">']
    mom = moments or []
    out.append(f'<UserMoments nb="{len(mom)}">' + "".join(
        f'<UserEntry id="{k}" content="{m}" window="{time_window}" />' for k, m in enumerate(mom)) + "</UserMoments>")
    out.append("</TimeSettings>")
    P = params or {}
    out.append(f'<Parameters nb="{len(P)}">')
    for k, (name, pts) in enumerate(P.items()):
        out.append(f'<Parameter id="{k}" name="{name}" logXinterp="0" logYinterp="0">' + "".join(
            f'<Moment id="{j}" t="{t:.9g}" value="{v:.9g}" />' for j, (t, v) in enumerate(pts)) + "</Parameter>")
    out += ["</Parameters>", "</MolflowSimuSettings>", "</SimulationEnvironment>"]
    pathlib.Path(path).write_text("\n".join(out))
    return dict(n_facets=n, n_vertices=len(V))


def run_cli(xml, ndes=None, time_s=None, outdir="molflow_res", threads=None, params_file=None, set_params=None,
            cli=None, log=None):
    """Ejecuta molflowCLI y devuelve la ruta del XML de resultados (con la tabla facet_physics.csv al lado)."""
    cli = cli or find_cli()
    if not cli:
        raise FileNotFoundError("molflowCLI no encontrado (define MOLFLOW_CLI o copia Molflow en <proyecto>/molflow/)")
    outdir = pathlib.Path(outdir); shutil.rmtree(outdir, ignore_errors=True)
    cmd = [cli, "-f", str(xml), "--outputPath", str(outdir), "--writeFacetQuantities", "-s", "0", "--noProgress",
           "-a", "0"]
    if ndes:
        cmd += ["-d", str(int(ndes))]
    if time_s:
        cmd += ["-t", str(int(time_s))]
    if threads:
        cmd += ["-j", str(int(threads))]
    if params_file:
        cmd += ["--setParamsByFile", str(params_file)]
    if set_params:
        cmd += ["--setParams", *set_params]
    r = subprocess.run(cmd, cwd=pathlib.Path(cli).parent, capture_output=True, text=True)
    if log:
        pathlib.Path(log).write_text(r.stdout + r.stderr)
    od = outdir if outdir.is_absolute() else pathlib.Path(cli).parent / outdir
    cands = sorted(od.glob("out_*")) if od.exists() else []
    if r.returncode != 0 or not cands:
        raise RuntimeError(f"molflowCLI falló ({r.returncode}):\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    res = cands[0]
    if res.suffix == ".zip":                               # entrada .zip → salida .zip: se extrae el XML
        import zipfile
        with zipfile.ZipFile(res) as z:
            name = next(n for n in z.namelist() if n.endswith(".xml"))
            z.extract(name, od)
        res = od / name
    return res


def read_results(out_xml, n_molecules=None, n_tp=None):
    """Resultados por momento: hits, absorbidos, desorbidos y presión por faceta (mbar).
    Flujo constante: escala de Molflow (molecules_per_TP en moléculas/s por partícula de prueba).
    Momentos (modo dependiente del tiempo): el atributo molecules_per_TP del XML NO está en moléculas por
    partícula; se usa n_molecules = moléculas reales inyectadas (∫Q dt / kT) y n_tp = partículas de prueba
    desorbidas (−d de molflowCLI). Verificado contra el analítico V/S y el TPMC (validation/molflow_vs_tpmc.py)."""
    root = ET.parse(out_xml).getroot()
    facets = root.find("Geometry").find("Facets")
    area = np.array([float(f.find("Area").get("m2")) for f in facets])
    mass = float(root.find("MolflowSimuSettings").find("Gas").get("mass"))
    out = dict(area_m2=area, mass=mass, moments=[])
    for m in root.find("MolflowResults").find("Moments").findall("Moment"):
        n = len(area)
        hit = np.zeros(n); absd = np.zeros(n); des = np.zeros(n); svo = np.zeros(n)
        for f in (m.find("FacetResults") if m.find("FacetResults") is not None else []):
            i = int(f.get("id")); h = f.find("Hits")
            hit[i] = float(h.get("nbHitEquiv")); absd[i] = float(h.get("nbAbsEquiv")); des[i] = float(h.get("nbDes"))
            svo[i] = float(h.get("sum_v_ort"))
        mptp = float(m.get("molecules_per_TP"))
        tw = float(m.get("timeWindow")) or 1.0
        if m.get("time") != "Constant flow" and n_molecules is not None:
            mptp = n_molecules / n_tp
        gl = m.find("Global"); g = gl.find("Hits") if gl is not None else None
        P = svo * mptp * mass * AMU_MOLFLOW / area / tw / 100.0              # Pa → mbar
        out["moments"].append(dict(time=m.get("time"), window=tw, molecules_per_TP=mptp, hits=hit, abs=absd,
                                   des=des, P_mbar=P, totalDes=float(g.get("totalDes")) if g is not None else float(des.sum())))
    return out
