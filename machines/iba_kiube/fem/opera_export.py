"""Paquete de modelo Opera-3d (Modeller + TOSCA magnetostático) del imán del KIUBE.

Genera out/opera3d/ con:
  geometria/            STEP por dominio (mm): sector de 90° (periódico) y mitad de 360°
  bh_acero_bajoC_kiube.csv  curva B-H calibrada (B [T], H [A/m]) para el editor de materiales de Opera
  bobina_solenoide.json valores del diálogo SOLENOID (XP1..YP4, CU1..CU4, densidad de corriente)
  OPERA_SETUP.md        pasos en el Modeller: importar, materiales, conductor, simetrías, malla, TOSCA
  referencia_ngsolve.json  resultado del mismo modelo resuelto aquí (NGSolve), para comparar
Y fem/opera_compare.py lee la tabla de Bz que exportes de Opera y calcula ⟨B⟩(r), flutter, B_colina y
B_valle frente a NGSolve, S09 y KEY22.

Las órdenes .comi de Opera no se escriben aquí: la sintaxis del Modeller cambia entre versiones y el
manual público disponible (v10, 2004) no la detalla; los valores de cada diálogo sí van completos.
"""
import json, pathlib, shutil, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import yaml
from kiube import load_params


def main(NI_opera=None):
    P, raw = load_params()
    M = yaml.safe_load(open(ROOT / "materials.yaml"))
    out = ROOT / "out/opera3d"; shutil.rmtree(out, ignore_errors=True); (out / "geometria").mkdir(parents=True)
    for tag in ("sector90", "half360"):
        src = ROOT / "out/fem_mag" / tag
        shutil.copytree(src / "domains", out / "geometria" / tag)
        shutil.copy(src / "domains.json", out / "geometria" / f"{tag}_domains.json")
    bh = M["bh_acero_bajoC_kiube"]
    np.savetxt(out / "bh_acero_bajoC_kiube.csv", np.c_[bh["B"], bh["H"]], delimiter=",", fmt="%.6g",
               header="B_T,H_A_m  (acero bajo C calibrado v0.18, J_s 2.12 T; extrapolar con mu0 por encima del ultimo punto)",
               comments="")
    b0 = M["bh_steel_1010"]
    np.savetxt(out / "bh_steel_1010_generica.csv", np.c_[b0["B"], b0["H"]], delimiter=",", fmt="%.6g",
               header="B_T,H_A_m  (1010 generica J_s 2.0 T, solo para comparar: pide 170 kA.v)", comments="")
    CO = P["coil"]
    ri, w, h, zb = CO["inner_radius"], CO["radial_width"], CO["height"], CO["z_bottom"]
    ref = ROOT / "out/fem_mag/ngsolve/resumen.json"
    NI_cal = None
    if ref.exists():
        R = json.load(open(ref))
        ni = np.array([r["NI_total"] for r in R]); b = np.array([r["Bavg_045"] for r in R])
        NI_cal = float(np.interp(1.3317, b, ni)) if len(R) > 1 else None   # ⟨B⟩(0.45 m) de KEY22 (13 317 G en el punto)
    NI = NI_opera or NI_cal or CO.get("NI_FE_calibrado", CO["NI_total"])
    J = NI / 2 / (w * h)                                                # A/mm² (por bobina)
    sol = dict(
        tipo="SOLENOID (Modeller > Create > Conductor > Solenoid)",
        sistema_local="origen (0, 0, 0), eje local Z = eje de la máquina; XP = radio, YP = z (mm)",
        bobina_superior=dict(XP1=ri, YP1=zb, XP2=ri, YP2=zb + h, XP3=ri + w, YP3=zb + h, XP4=ri + w, YP4=zb,
                             CU1=0, CU2=0, CU3=0, CU4=0),
        bobina_inferior="espejo en z = 0 (o simetría de plano medio en el modelo de mitad)",
        densidad_corriente_A_mm2=round(J, 4), NI_total_A_vuelta=NI,
        NI_origen=("calibrado con NGSolve para ⟨B⟩ = 13 253 G en 0.38–0.47 m (KEY22)" if NI_cal and not NI_opera
                   else "params.yaml: NI_FE_calibrado (v0.18; 105 kA·v por potencia, IPAC26)"),
        nota="signo de J: invertir si Bz sale negativo; H⁻ gira en sentido antihorario con Bz > 0")
    json.dump(sol, open(out / "bobina_solenoide.json", "w"), indent=1, ensure_ascii=False)
    if (ROOT / "out/fem_mag/ngsolve").exists():
        shutil.copytree(ROOT / "out/fem_mag/ngsolve", out / "referencia_ngsolve")
    (out / "OPERA_SETUP.md").write_text(f"""# Cyclone KIUBE — modelo Opera-3d (Modeller + TOSCA magnetostático)

El mismo modelo se resolvió aquí con NGSolve (`fem/magnet_ngsolve.py`, misma física): los resultados
están en `referencia_ngsolve/` para comparar el tuyo con `python fem/opera_compare.py tu_tabla.csv`.

## 1. Geometría
- `File > Import` → `geometria/sector90/*.step` (mm). Son dominios conformes: hierro (yugo A+B, sectores C,
  polos D, insertos E, tapón F), bobina y aire. En Opera el aire suele ser la *background region*: puedes omitir
  `air.step` y crear el fondo con *Model > Background* (≈1.6 × el lado del yugo, como el STEP).
- Mitad superior z ≥ 0; sector de 90° entre los centros de colina (45° y 135°).
- Geometría v0.18: tapa A de 181 mm, valle con piso a 95 mm del plano medio y tapón rebajado 16 mm (ajustados con este mismo FE contra KEY22).

## 2. Materiales
- Hierro: material no lineal con la curva `bh_acero_bajoC_kiube.csv` (B [T], H [A/m]): acero de bajo C con J_s = 2.12 T,
  calibrada en v0.18 (CALC). Reemplázala por la curva del lote si la consigues. `bh_steel_1010_generica.csv` es solo para comparar.
- Bobina y aire: μr = 1. El cobre de dees, el inox de la cámara y el Nb de los blancos no son magnéticos: se omiten.

## 3. Conductor (bobina)
`Create > Conductor > Solenoid` con los valores de `bobina_solenoide.json`:
- XP1 = {ri:.0f}, YP1 = {zb:.0f}; XP2 = {ri:.0f}, YP2 = {zb + h:.0f}; XP3 = {ri + w:.0f}, YP3 = {zb + h:.0f}; XP4 = {ri + w:.0f}, YP4 = {zb:.0f}; CU1–CU4 = 0
- Densidad de corriente: **{J:.4f} A/mm²** (NI total = {NI:,.0f} A·v, {sol['NI_origen']}).

> **NI y curva B-H (v0.18).** Con la curva 1010 genérica (J_s 2.0 T) el FE necesitaba ≈ 170 kA·v. Con el acero de bajo C
> calibrado (J_s 2.12 T), la tapa de 181 mm y el valle a 95 mm necesita ≈ 116 kA·v, frente a ≈ 105 kA·v de la potencia (IPAC26).
> Si tienes la curva del lote, cárgala y recalibra J hasta ⟨B⟩ = 13 253 G en 0.38–0.47 m; el flutter debe quedar en ≈ 0.39.

## 4. Condiciones de borde y simetrías
- z = 0: **Normal magnetic** (campo normal al plano medio; simetría par del dipolo).
- φ = 45° y φ = 135°: **simetría rotacional periódica** de 90° (o modela la mitad de 360° de `geometria/half360/`).
- Borde exterior del fondo: **Tangential magnetic** (flujo tangente).

## 5. Malla (TOSCA)
- Caras de polo, inserto y tapón: tamaño máximo 3–5 mm; el gap de colina de 24 mm debe tener ≥ 6 elementos en z.
- Bordes del escalón 2 y cortes correctores: ≤ 2 mm.
- Sectores C: 10–15 mm; yugo: 40–60 mm; fondo: crecimiento 1.3–1.5, hasta 250 mm.
- Elementos cuadráticos recomendados. Referencia de convergencia aquí: con 7 mm en el gap (lineales) ⟨B⟩ cambia
  < 0.3 % frente a 10 mm.

## 6. Solución y exportación
- TOSCA magnetostático no lineal (Newton), tolerancia 1e-4.
- Exporta Bz en el plano medio (z = 0.2 mm) en una grilla polar r = 0…600 mm (paso 5 mm), φ = 45°…135° (paso 0.5°),
  como tabla de texto/CSV con columnas X, Y, Z, BZ (mm, T o gauss).
- Compara: `python fem/opera_compare.py tabla_opera.csv` → ⟨B⟩(r), flutter, B_colina, B_valle frente a NGSolve,
  SUPAG09 (2.15 / 0.45 T) y KEY22 (⟨B⟩ y F del KIUBE).
""")
    print("ok", out, "NI", NI)


def region_central(preset="c18_9"):
    """Paquete Opera-3d (Modeller + TOSCA electrostático) de la región central desde fem/central_region.py."""
    src = ROOT / f"out/central_region/{preset}"
    R = json.load(open(src / "resultado.json")); pr = R["parametros"]
    out = ROOT / f"out/opera3d/region_central_{preset}"; shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True)
    if (src / "electrodos").exists():
        shutil.copytree(src / "electrodos", out / "electrodos")
    if (src / "campo_E_referencia_FE.txt").exists():
        shutil.copy(src / "campo_E_referencia_FE.txt", out / "campo_E_referencia_FE.txt")
    fu = "\n".join(f"- {F['especie']} (h = {F['h']}): chimenea en ({F['S'][0]:.2f}, {F['S'][1]:.2f}) mm, rendija hacia "
                   f"({F['e_n'][0]:.3f}, {F['e_n'][1]:.3f})" for F in R["fuentes"])
    (out / "OPERA_REGION_CENTRAL.md").write_text(f"""# Región central {preset} — modelo Opera-3d (TOSCA electrostático)

Mismo modelo que `fem/central_region.py` (NGSolve, H1 de orden 2, {R['malla']['ne']:,} tetraedros, {R['malla']['ndof']:,} gdl).
Parametrizado: cambia cualquier valor con `python fem/central_region.py --preset {preset} --set clave=valor` y se
regeneran los STEP de `electrodos/`.

## 1. Geometría (`electrodos/*.step`, mm)
- `dees_*`, `pullers_*`: electrodo de RF. `colinas_tapon_*`, `dummy_dees_*`, `chimeneas_*`: tierra.
- Fondo (vacío): cilindro r = {pr['R_dom']:.0f} mm, |z| ≤ {pr['H_dom']:.0f} mm (estructura completa de los dees).
- Fuentes:
{fu}

## 2. Potenciales y contorno
- Dees + pullers: **1 V** (el tracking escala por V_dee·cos(ωt + φ₀) con V_dee = {pr['V_dee'] / 1e3:.0f} kV; cuasi-estático, λ_RF ≈ 7 m).
- Colinas, tapón, dummy dees, chimeneas y borde del fondo: **0 V**.
- z = 0 puede usarse como plano de simetría (**Tangential electric**: E_z = 0) para modelar la mitad.

## 3. Malla (graduada, como en el modelo de referencia)
| Zona | Tamaño |
|---|---|
| Gap fuente–puller y rendija | {pr['h_gap']} mm |
| Narices de los dees y centro (r < {pr.get('R_fino', 40):.0f} mm) | {pr['h_nose']} mm |
| r < 150 mm | {pr.get('h_mid', pr['h_far'])} mm |
| Resto de la estructura de los dees | {pr['h_far']} mm |
Elementos cuadráticos.

## 4. Exportación y uso de tu campo
- Exporta E en la grilla de `campo_E_referencia_FE.txt` (x, y = −60…60 mm cada 1 mm; z = −6…6 cada 1.5 mm),
  columnas `X Y Z EX EY EZ` en mm y V/mm (por 1 V de dee).
- Compara con la referencia FE y rehaz las órbitas con TU campo:
  `python fem/central_region.py --preset {preset} --campoE tu_tabla_opera.txt`
  (el tracking relativista, el campo magnético y las figuras son los mismos).
""")
    print("ok", out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "region_central":
        region_central(sys.argv[2] if len(sys.argv) > 2 else "c18_9")
    else:
        main()
