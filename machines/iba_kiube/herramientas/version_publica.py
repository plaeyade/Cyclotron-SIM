"""Genera la VERSIÓN PÚBLICA del proyecto KIUBE, sin datos de documentos propietarios o confidenciales de IBA
(Site Requirement Guide 'SRG' y cotización del proyecto 'PROP').

Qué hace:
  1. Copia el proyecto (sin out/, cachés ni binarios) a dist/kiube_publico/.
  2. params.yaml: cada cota SRG/PROP vuelve a su valor público anterior (DS, DSR, TUD03, EST) o se elimina si no
     tenía equivalente público; se eliminan los bloques 'proyecto', 'sitio' y 'fuente_rad'.
  3. materials.yaml: notas sin referencias al SRG; el autoblindaje usa una composición genérica (EST).
  4. Excluye los archivos derivados del SRG (benchmark de blindaje, validar_srg, registro de fuentes completo) y limpia
     README y comentarios.
  5. Comprueba: no quedan referencias a SRG/PROP/cliente, y el modelo construye (pytest rápido en la copia).
Uso: python herramientas/version_publica.py [--sin_pruebas]"""
import pathlib, re, shutil, subprocess, sys
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
DST = ROOT / "dist/kiube_publico"
EXCLUIR_DIRS = {"out", "dist", "__pycache__", "mc_benchmark", ".pytest_cache", "build"}
EXCLUIR_ARCH = {"validation/validar_srg.py", "validation/validacion_srg.json", "validation/benchmark_srg.py",
                "validation/benchmark_srg.json", "validation/benchmark_srg.png", "docs/FUENTES.md",
                "herramientas/version_publica.py"}
# (grupo, cota) → valor público {v, src, note} o None (eliminar)
PUBLICO = {
    ("yoke", "floor_to_top"): dict(v=1550.0, src="JUAS22", note="plano medio a 1.12 m del piso (JUAS 2022, tabla 18P) + H/2"),
    ("ion_source", "actuator_len"): dict(v=170.0, src="DSR"),
    ("extraction", "strippers_instalados"): None,
    ("target", "types"): dict(v=["F18"] * 8, src="EST", note="F18 | N13 | C11 | Ga68 | none"),
    ("target", "volumen_ml"): dict(v={"F18": 3.6, "Ga68": 3.7, "C11": 55.0}, src="CBR",
                                   note="F18 ficha DS; Ga68 y 11C folleto público CBR"),
    ("target", "puertas_blindadas"): None,
    ("vault", "ceiling"): dict(v=1500.0, src="EST"),
    ("exterior", "overall_xyz"): dict(v=[1950.0, 1950.0, 1850.0], src="JUAS22", note="JUAS 2022, tabla 18P (ficha DS: 1.9 × 1.9 × 1.8)"),
    ("exterior", "weight_t"): dict(v=17.0, src="JUAS22", note="JUAS 2022 (ficha DS: 18 t)"),
    ("exterior", "leg_square"): dict(v=200.0, src="DSR"),
    ("exterior", "leg_plate"): None,
    ("exterior", "leg_xy"): dict(v=560.0, src="DSR"),
    ("exterior", "eyebolt"): None,
    ("exterior", "self_shield"): dict(v=False, src="DS"),
    ("exterior", "self_shield_piece"): dict(v=[4645.0, 2770.0, 2575.0], src="DS"),
    ("exterior", "self_shield_mass_kg"): None,
    ("exterior", "self_shield_wall"): dict(v=410.0, src="EST"),
}
BLOQUES_FUERA = {"proyecto", "sitio", "fuente_rad"}
PATRON = re.compile(r"\bSRG\b|\bPROP\b|Zumaya|55672|2025/8499|Site Requirement")
# reemplazos de comentarios/textos en código (no cambian el comportamiento)
TEXTOS = {
    "kiube/mc.py": [("    # PROP: KIUBE 180", "    # bloque 'proyecto' si existe"),
                    ("'SRG' (IBA 6.9e-3 n/p, mismo espectro)", "'SRG' (fuente_rad.n_por_p_18O, si existe)")],
    "kiube/geometry.py": [("   # PROP: 4 de serie", "   # por defecto, los 8"),
                          ("# pie 350 × 350 (SRG-0315)", "# placa de apoyo (exterior.leg_plate)"),
                          ("# ---- cáncamo central M48 (SRG): único punto de izado superior; su cima marca la altura total (1850 mm)",
                           "# ---- cáncamo central de izado (opcional: exterior.eyebolt)"),
                          ("# CALC: 48.1 t por mitad (SRG) con ρ ≈ 3.2 g/cm³", "# EST")],
    "tests/test_modelo.py": [("# ------------------------------------------------------------------ documentos (SRG) — solo si los datos están",
                              "# ------------------------------------------------------------------ documentos de instalación — solo si existen"),
                             ('pytest.skip("versión pública: sin datos del SRG")', 'pytest.skip("versión pública: sin datos de instalación")')],
}


def copiar():
    if DST.exists():
        shutil.rmtree(DST)
    for p in ROOT.rglob("*"):
        rel = p.relative_to(ROOT)
        if any(part in EXCLUIR_DIRS for part in rel.parts) or str(rel) in EXCLUIR_ARCH or p.is_dir():
            continue
        if p.suffix in (".step", ".h5m", ".gdml", ".zip") or p.name == "srg_bench":
            continue
        q = DST / rel; q.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, q)


def params_publicos():
    """Reescribe params.yaml línea a línea (conserva comentarios y orden)."""
    src = (ROOT / "params.yaml").read_text().splitlines()
    out, grupo, saltar = [], None, False
    for ln in src:
        mg = re.match(r"^([a-z_0-9]+):", ln)
        if mg:
            grupo = mg.group(1); saltar = grupo in BLOQUES_FUERA
            if saltar:
                continue
        if saltar:
            if ln.strip() == "" or ln.startswith(" ") or ln.startswith("#"):
                continue
        mk = re.match(r"^  ([a-z_0-9]+):\s+\{", ln)
        if mk and (grupo, mk.group(1)) in PUBLICO:
            nv = PUBLICO[(grupo, mk.group(1))]
            if nv is None:
                continue
            item = yaml.safe_dump({mk.group(1): nv}, default_flow_style=True, allow_unicode=True, width=1000).strip()
            item = item[1:-1]                                  # quita las llaves exteriores
            ln = "  " + item
        if ln.lstrip().startswith("#") and PATRON.search(ln):  # cabecera: quita las etiquetas SRG/PROP
            continue
        if PATRON.search(ln):
            ln = re.sub(r",?\s*note:\s*\"[^\"]*\"", "", ln) if mk else ln
        out.append(ln)
    txt = "\n".join(out) + "\n"
    txt = txt.replace("(v0.6)", "(versión pública)")
    (DST / "params.yaml").write_text(txt)
    yaml.safe_load(txt)                                       # valida que sigue siendo YAML


def materiales_publicos():
    m = (DST / "materials.yaml").read_text()
    m = re.sub(r"ρ = 2\.35 = mínimo exigido por IBA \(SRG §2\.1\.2; antes 2\.30\)", "ρ = 2.35 (hormigón de blindaje habitual)", m)
    m = re.sub(r'(selfshield_magnetita:\n  src: )"[^"]*"\n  density: [0-9.]+\n  wt: \{[^}]*\}',
               r'\1"autoblindaje: mezcla genérica magnetita + agua (EST)"\n  density: 3.2\n'
               r'  wt: {Fe: 62.0, O: 35.5, H: 2.0, Si: 0.3, Ca: 0.2}', m)
    (DST / "materials.yaml").write_text(m)


def textos():
    for f, pares in TEXTOS.items():
        p = DST / f
        t = p.read_text()
        for a, b in pares:
            t = t.replace(a, b)
        p.write_text(t)
    # README: sin aviso confidencial ni secciones v0.12/v0.13; líneas con SRG/PROP fuera
    r = (DST / "README.md").read_text()
    r = re.sub(r"\n> \*\*Desde v0\.12 incluye datos[^\n]*\n", "\n", r)
    r = re.sub(r"## Qué cambió en v0\.1[23]\n.*?(?=\n## )", "", r, flags=re.S)
    r = "\n".join(l for l in r.splitlines() if not PATRON.search(l)) + "\n"
    r = r.replace("Modelo reconstruido a partir de literatura pública y de los documentos aportados.",
                  "Modelo reconstruido solo a partir de literatura pública. **Versión pública**: sin datos de documentos de IBA.")
    (DST / "README.md").write_text(r)
    rf = DST / "fem/README_RF.md"
    rf.write_text("\n".join(l for l in rf.read_text().splitlines() if not PATRON.search(l)) + "\n")
    (DST / "docs/FUENTES.md").write_text(
        "# Fuentes\n\nVersión pública: todas las cotas salen de literatura pública (incluida la clase JUAS 2022 de IBA), de la "
        "ficha técnica y el folleto de IBA (públicos) o son estimaciones (EST) y cálculos (CALC). Ver la cabecera de "
        "`params.yaml`.\n")


def auditar():
    restos = []
    for p in DST.rglob("*"):
        if p.is_file() and p.suffix in (".py", ".md", ".yaml", ".json", ".txt", ".sh", ".csv"):
            for i, l in enumerate(p.read_text(errors="ignore").splitlines(), 1):
                if PATRON.search(l) and not (p.name == "test_modelo.py" and "ETIQUETAS_EXTRA" in l):
                    restos.append(f"{p.relative_to(DST)}:{i}: {l.strip()[:100]}")
    return restos


if __name__ == "__main__":
    copiar(); params_publicos(); materiales_publicos(); textos()
    restos = auditar()
    print(f"[publico] copia en {DST}")
    print(f"[publico] referencias restantes: {len(restos)}")
    for r in restos:
        print("   ", r)
    if "--sin_pruebas" not in sys.argv:
        r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=DST,
                           capture_output=True, text=True)
        print("[publico] pruebas:", r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-500:])
        sys.exit(1 if restos or r.returncode else 0)
    sys.exit(1 if restos else 0)
