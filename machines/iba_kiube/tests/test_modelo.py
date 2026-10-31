"""Batería de pruebas del modelo KIUBE.   Rápidas:  python -m pytest      Todas:  python -m pytest -m ""
Cada prueba falla si un cambio rompe una propiedad verificada antes (fuentes, geometría, masas, física)."""
import json, math, re, subprocess, sys
import numpy as np
import pytest
from conftest import ROOT

ETIQUETAS_EXTRA = {"EST", "CALC", "DS", "DSR", "F7", "RPOS", "SRG", "CBR", "PROP"}


# ------------------------------------------------------------------ parámetros y materiales
def test_toda_cota_tiene_fuente_declarada(PE):
    _, _, raw = PE
    cab = open(ROOT / "params.yaml").read().split("\n\n")[0]
    tags = set(re.findall(r"^#\s+([A-Z][A-Z0-9-]+)\s", cab, re.M)) | ETIQUETAS_EXTRA
    sin = [(g, k, v.get("src")) for g, d in raw.items() for k, v in d.items() if v.get("src") not in tags]
    assert not sin, f"cotas con fuente no declarada en la cabecera: {sin}"


def test_materiales_definidos(lods, mats):
    usados = {p["mat"] for ps in lods.values() for p in ps}
    falt = [m for m in usados if m not in mats or "density" not in mats[m]]
    assert not falt, falt


# ------------------------------------------------------------------ geometría
@pytest.mark.parametrize("lod", ["print", "neutronics", "fem_struct"])
def test_geometria_valida_sin_solapes(lods, lod):
    import build as B
    ps = lods[lod]
    assert len({p["name"] for p in ps}) == len(ps), "nombres de pieza repetidos"
    inval = [p["name"] for p in ps if not p["shape"].isValid()]
    assert not inval, inval
    assert B.overlaps(ps) == []


def test_masa_frente_al_srg(lods, mats, PE):
    P, _, _ = PE
    m = sum(p["shape"].Volume() * 1e-6 * mats[p["mat"]]["density"] for p in lods["fem_struct"])
    ref = P["exterior"]["weight_t"] * 1000
    assert abs(m / ref - 1) < 0.10, f"masa {m:.0f} kg frente a {ref:.0f} kg"


def test_cotas_de_instalacion(lods, PE):
    P, _, _ = PE
    H = P["yoke"]["height"] / 2
    floor = H - P["yoke"]["floor_to_top"]
    zmax = max(p["shape"].BoundingBox().zmax for p in lods["print"])
    zmin = min(p["shape"].BoundingBox().zmin for p in lods["print"])
    assert abs(zmin - floor) < 1.0, "la máquina no apoya en el piso"
    tol = 40.0 if P["exterior"].get("eyebolt") else 80.0   # sin el cáncamo central (versión pública) faltan ~70 mm
    assert abs((zmax - floor) - P["exterior"]["overall_xyz"][2]) < tol, f"altura total {zmax - floor:.0f} mm"


def test_extraccion_8_puertos(PE):
    _, E, _ = PE
    assert len(E["targets"]) == 8 and len(E["strippers"]) == 8


# ------------------------------------------------------------------ vacío (TPMC propio) frente a Clausing
@pytest.mark.parametrize("LR", [1.0, 10.0])
def test_tpmc_transmision_clausing(LR):
    from kiube.tpmc import TPMC, tube_mesh, CLAUSING
    m, c1, c2 = tube_mesh(R=1.0, L=LR, n=64)
    s = np.where(c1 | c2, 1.0, 0.0)
    N = 40000
    r = TPMC(m, s, seed=7).run(np.where(c1, m.area_faces, 0.0), n_particles=N)
    W = r["absorbed"][c2].sum() / r["n"]
    sig = math.sqrt(CLAUSING[LR] * (1 - CLAUSING[LR]) / N)
    assert abs(W - CLAUSING[LR]) < 4 * sig + 0.003, f"W = {W:.4f} vs {CLAUSING[LR]}"


# ------------------------------------------------------------------ paleta y salidas de impresión
def test_paleta_cubre_todas_las_piezas(lods):
    from kiube import palette as PAL
    mats = PAL.cargar()["materiales"]
    for p in lods["print"]:
        for v in ("foto", "render"):
            assert PAL.material_visual(p["name"], v) in mats, p["name"]


def test_stl_hd_estanco_si_existe():
    man = ROOT / "out/stl_hd/escala_1-10_alta/manifest.json"
    if not man.exists():
        pytest.skip("sin salida stl_hd (python build.py --lod stl_hd)")
    m = json.load(open(man))
    r = m["resumen"]
    assert r["estancos"] == r["cuerpos"] == r["manifold"] and r["con_autointersecciones"] == 0
    assert m["completo"]["estanca"] and m["completo"]["componentes"] == 1


def test_3mf_valido_si_existe():
    f = ROOT / "out/stl_hd/escala_1-10_alta/completo/kiube_completo.3mf"
    if not f.exists():
        pytest.skip("sin 3MF")
    lib3mf = pytest.importorskip("lib3mf")
    w = lib3mf.get_wrapper(); m = w.CreateModel(); rd = m.QueryReader("3mf"); rd.SetStrictModeActive(True)
    rd.ReadFromFile(str(f))
    assert rd.GetWarningCount() == 0
    it = m.GetMeshObjects(); bad = 0
    while it.MoveNext():
        bad += 0 if it.GetCurrentMeshObject().IsManifoldAndOriented() else 1
    assert bad == 0


# ------------------------------------------------------------------ documentos (SRG) — solo si los datos están
def test_validacion_srg(PE):
    P, _, _ = PE
    if "sitio" not in P:
        pytest.skip("versión pública: sin datos del SRG")
    r = subprocess.run([sys.executable, str(ROOT / "validation/validar_srg.py")], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-1500:]
    res = json.load(open(ROOT / "validation/validacion_srg.json"))["resumen"]
    assert res["ok"] == res["total"], res




# ------------------------------------------------------------------ imán (v0.18): calibración FE frente a KEY22/S09
def test_calibracion_iman(PE):
    """El caso FE calibrado (out/fem_mag/ngsolve, NI_FE_calibrado) reproduce KEY22 y S09, y el NI queda cerca de la potencia."""
    import yaml
    P, _, _ = PE
    f = ROOT / "out/fem_mag/ngsolve/resumen.json"
    if not f.exists():
        pytest.skip("sin resultado FE del imán (python fem/magnet_ngsolve.py)")
    sys.path.insert(0, str(ROOT / "validation"))
    import calibrar_iman as ci
    c = ci.caso(next((ROOT / "out/fem_mag/ngsolve").glob("resultado_NI*.json")))
    assert abs(c["NI_total"] - P["coil"]["NI_FE_calibrado"]) < 1
    # ⟨B⟩ 0.38–0.47 m en el nivel de operación isócrono de cada mapa (v0.34: el FE con shim se opera a ×escala_iso,
    # el campo medido de KEY22 a ×1.0038 con f_RF 40.640 MHz)
    k = P["coil"].get("escala_iso", 1.0)
    assert abs(c["B_banda_T"] * k / (ci.KEY["B_banda_T"] * (1.0038 if k != 1.0 else 1.0)) - 1) < 0.005
    assert abs(c["F_banda"] - ci.KEY["F_banda"]) < 0.03                    # flutter
    assert abs(c["B_centro_rel"] - ci.KEY["B_centro_rel"]) < 0.03          # joroba central
    # S09. v0.34: el escalón 2 a 526 mm y el shim bajan el valle a 0.42 T en el nivel de operación (−7 % frente a S09)
    assert abs(c["B_valle_045"] * k - 0.45) < 0.04 and abs(c["B_colina_045"] * k - 2.15) < 0.10
    assert P["coil"]["NI_FE_calibrado"] / P["coil"]["NI_total"] < 1.15     # ≤ 15 % sobre la potencia
    M = yaml.safe_load(open(ROOT / "materials.yaml"))
    assert M["fem_properties"]["steel_1010"]["magnetic"] == "bh_acero_bajoC_kiube"
    bh = M["bh_acero_bajoC_kiube"]; J = np.array(bh["B"]) - 4e-7 * math.pi * np.array(bh["H"])
    assert 2.0 < J[-1] <= 2.16                                            # J_s físico (hierro puro 2.16 T)


# ------------------------------------------------------------------ vacío (v0.20): bombeo frente al medido de IBA
def test_bombeo_frente_a_tud03(PE):
    """Con el acero calibrado, el pump-down a 30 y 45 min queda a ±50 % del medido (TUD03T, 4 × VHS-6)."""
    P, _, _ = PE
    f = ROOT / "out/vacuum/tpmc_pumpdown.json"
    pm = P["vacuum"].get("pumpdown_medido")
    if not f.exists() or not pm:
        pytest.skip("sin resultado TPMC del vacío o sin datos medidos")
    from kiube.pumpdown import curve
    T = json.load(open(f))
    _, Pm = curve(P, T, t_h=np.array(pm["t_min"][1:]) / 60)
    for mod, med in zip(Pm, pm["P_mbar"][1:]):
        assert 0.5 < mod / med < 1.5, (mod, med)
