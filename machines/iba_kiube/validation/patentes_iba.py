"""Cotas del modelo frente a los rangos de las 4 patentes de IBA de 2016 (prioridad EP16169489/-90/-94/-97, 13-05-2016):
US9907153 (ciclotrón compacto = «geometry»), US10064264 (inserto de polo), US9961757 (borde de colina = «pole contour»),
US10278277 (corrector de gradiente). Ninguna da cotas absolutas del KIUBE ni potencia/NI: solo rangos y reivindicaciones.
Uso: python validation/patentes_iba.py → validation/patentes_iba.json"""
import json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from kiube import load_params

P, _ = load_params()
Y, PO, IS = P["yoke"], P["pole"], P["ion_source"]
O = json.load(open(ROOT / "out/neutronics/report.json"))["optics"]
Gh = PO["gap"]; Gv = 2 * PO["valley_floor_z"]
Ra, La = Y["pump_hole_radius"], Y["pump_hole_r_pos"]
Lv = PO["sector_C_radius"]                                  # borde periférico del suelo del valle
Tv = Y["side"] / 2 - Y["cavity_radius"]                     # retorno de flujo en el centro de la cara (valle)
Lh = PO["step2_edge_L"]; ah = math.radians(PO["hill_angle"]); Ah = Lh * ah
av = 360 / PO["n_sectors"] - PO["hill_angle"]
Lap = Lv - (La + Ra)
Hh = PO["valley_floor_z"] - Gh / 2                          # altura de la colina sobre el suelo del valle

def chk(nombre, valor, lo, hi, patente, tipo="preferido"):
    return dict(cota=nombre, modelo=round(valor, 4), rango=[lo, hi], patente=patente, tipo=tipo,
                cumple=(lo is None or valor >= lo) and (hi is None or valor <= hi))

C = [
    chk("Gh gap de colina [mm]", Gh, 22, 26, "US9907153"),
    chk("Gv gap de valle [mm]", Gv, 200, 250, "US9907153 reiv. 19", "reivindicación"),
    chk("Gh/Gv", Gh / Gv, 0.08, 0.20, "US9907153"),
    chk("αv ángulo de valle [°]", av, 42, 45, "US9907153 reiv. 20", "reivindicación"),
    chk("αh ángulo de colina [°]", PO["hill_angle"], 44, 50, "US9907153"),
    chk("Gh/tan(αv) [mm]", Gh / math.tan(math.radians(av)), None, 30, "US9907153"),
    chk("2Ra diámetro del agujero de bombeo [mm]", 2 * Ra, 240, 300, "US9907153"),
    chk("2Ra/Lv", 2 * Ra / Lv, 0.48, 0.55, "US9907153"),
    chk("2Ra/La", 2 * Ra / La, 0.70, None, "US9907153 reiv. 16", "reivindicación"),
    chk("Lap = Lv − (La + Ra) [mm]", Lap, 1, 30, "US9907153"),
    chk("(Lap·Tv)/Lv²", Lap * Tv / Lv ** 2, None, 0.03, "US9907153"),
    chk("(Gv·Tv)/Lv²", Gv * Tv / Lv ** 2, None, 0.10, "US9907153 reiv. 17", "reivindicación"),
    chk("(Gh·Gv)/Lv²", Gh * Gv / Lv ** 2, None, 0.02, "US9907153 reiv. 18", "reivindicación"),
    chk("Lh radio del borde de colina [mm]", Lh, 520, 550, "US9907153 / US9961757 («típico»)"),
    chk("Rh/Lh", O["Rh_over_Lh"], None, 0.65, "US9961757"),
    chk("W9 ancho del inserto [mm]", PO["insert_groove_w"], 15, 150, "US10064264"),
    chk("W9/Ah", PO["insert_groove_w"] / Ah, None, 0.15, "US10064264"),
    chk("H9/Hh altura del inserto / altura de colina", PO["insert_groove_d"] / Hh, 0.40, 0.60, "US10064264"),
    chk("H10/Lh profundidad del corrector de gradiente", PO["gc_cut_depth"] / Lh, 0.06, 0.10, "US10278277"),
    chk("L10/Gh ancho del corrector / gap", PO["gc_cut_width"] / Gh, 3.0, 3.5, "US10278277"),
]
# relación del borde con la órbita (US9961757): Rh ≈ Rp + 0.6·Gh, con el centro del arco en el centro de curvatura de la órbita
rel = dict(edge_R_modelo=O["edge_R"], rho_orbita_colina=O["rho_hill"], Rh_patente=O["rho_hill"] + 0.6 * Gh,
           centro_arco_modelo=O["edge_d0"], centro_curvatura_orbita=O["r_ext_hill"] - O["rho_hill"])
out = dict(comprobaciones=C, borde_vs_orbita=rel,
           nota="Lv = radio del sector C (borde del suelo del valle); Tv = lado/2 − radio de la cavidad (cara, en el valle); "
                "Hh = suelo del valle − Gh/2. Las patentes no dan potencia del imán ni NI.")
json.dump(out, open(ROOT / "validation/patentes_iba.json", "w"), indent=1, ensure_ascii=False)
for c in C:
    r = f"[{c['rango'][0] if c['rango'][0] is not None else ''} – {c['rango'][1] if c['rango'][1] is not None else ''}]"
    print(f"{'OK ' if c['cumple'] else 'NO '} {c['cota']:<48} {c['modelo']:>9}  {r:<14} {c['tipo']}  {c['patente']}")
print(json.dumps(rel))
