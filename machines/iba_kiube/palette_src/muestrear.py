"""Genera palette.yaml (paleta de colores del KIUBE) muestreando las imágenes de las fichas técnicas IBA.

Método (reproducible; regiones en palette_src/regiones.py, figura de control en docs/paleta_regiones_*.png):
  0. Las imágenes CMYK de los PDF se pasan a sRGB con el perfil FOGRA39 (sin esto el cobre sale violeta).
  1. Mediana robusta de cada región (se descarta el 25 % más claro y el 25 % más oscuro: brillos y sombras).
  2. Corrección de 2 puntos en luz lineal por imagen: el yugo blanco y un hueco en sombra (negro) se suponen
     neutros. Elimina el viraje azul del etalonaje de las fotos (no lo elimina del todo en las sombras).
     Metales/colores saturados: se toma la mitad más cromática de la región (sin reflejos blancos).
  3. Color plano para imprimir/visualizar (sRGB) = muestra balanceada.
  4. Color base PBR (lineal) para Blender: pinturas y plásticos = muestra linealizada y escalada para que el
     blanco del yugo tenga albedo 0.80; metales = reflectancia F0 de tablas PBR estándar (el color de una
     foto de metal es el reflejo del entorno, no el albedo), con el tono verificado contra la foto.
  5. RAL más cercano por ΔE76 en CIELAB contra una tabla RAL Classic aproximada (sRGB de referencia): ORIENTATIVO,
     las fichas no indican RAL.
Uso: python palette_src/muestrear.py   (escribe palette.yaml y docs/paleta_kiube.png)"""
import json, pathlib, sys
import numpy as np
import yaml
from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "palette_src"))
from regiones import R, REF, CROMA                                           # noqa: E402

DOCS = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path("/home/claude/docs5")
FUENTE = {"v6_p1_img0.jpeg": "VE24 p.1 (render oficial)", "v6_p7_img0.jpeg": "VE24 p.7 (foto máquina real)",
          "k2_p4_img0.jpeg": "KEY22 p.4 (foto interior)"}
WHITE = {"v6_p1_img0.jpeg": "yugo_blanco_iluminado", "v6_p7_img0.jpeg": "yugo_blanco_foto"}

# RAL Classic: valores sRGB aproximados de referencia (varían según la fuente ±5/255)
RAL = {"9003 blanco señales": (236, 236, 231), "9010 blanco puro": (241, 236, 225), "9016 blanco tráfico": (241, 240, 234),
       "9002 blanco grisáceo": (215, 215, 205), "7035 gris luminoso": (203, 208, 204), "7047 telegris 4": (208, 208, 208),
       "7040 gris ventana": (157, 161, 170), "7001 gris plata": (138, 149, 151), "7031 gris azulado": (91, 104, 109),
       "7016 gris antracita": (56, 62, 66), "7024 gris grafito": (71, 74, 81), "7021 gris negruzco": (47, 50, 52),
       "9005 negro intenso": (14, 14, 16), "9006 aluminio blanco": (161, 161, 160), "9007 aluminio gris": (135, 133, 129),
       "5014 azul paloma": (96, 111, 140), "5024 azul pastel": (93, 155, 155), "5012 azul luminoso": (0, 137, 182),
       "5015 azul cielo": (34, 113, 179), "5005 azul señales": (30, 45, 110), "3020 rojo tráfico": (204, 6, 5),
       "3000 rojo vivo": (175, 43, 30), "1036 oro perlado": (112, 83, 53), "8023 anaranjado pardo": (166, 94, 46),
       "8004 pardo cobrizo": (143, 78, 53), "1011 beige pardo": (138, 102, 66), "7032 gris guijarro": (184, 183, 153),
       "7044 gris seda": (202, 196, 176)}

# PBR F0 (lineal) de metales: tablas estándar (Hoffman, «Physics and Math of Shading», SIGGRAPH; Real-Time Rendering)
F0 = {"hierro": (0.56, 0.57, 0.58), "cobre": (0.955, 0.638, 0.538), "aluminio": (0.913, 0.922, 0.924),
      "laton": (0.910, 0.778, 0.423), "inox": (0.60, 0.59, 0.57), "niobio": (0.58, 0.58, 0.60)}

# ---------------------------------------------------------------- especificación de materiales visuales
# foto/render: muestra 'imagen:region'; tipo pintura|metal|plastico|otro; pbr: clave F0 o None (→ de la muestra)
M = {
 "pintura_blanca":   dict(d="Pintura del yugo (caras exteriores), tapa de la fuente", tipo="pintura",
                          foto="v6_p7_img0.jpeg:yugo_blanco_foto", render="v6_p1_img0.jpeg:yugo_blanco_iluminado",
                          rough=0.30, coat=0.4, fil="PLA blanco (blanco frío)", conf="alta"),
 "banda_chaflan":    dict(d="Banda de los chaflanes de 45° del yugo", tipo="pintura",
                          foto="v6_p7_img0.jpeg:banda_chaflan_foto", render="v6_p1_img0.jpeg:banda_chaflan_render",
                          rough=0.35, fil="PLA gris azulado claro / gris plata", conf="media",
                          nota="Foto: pintura gris azulada. Render: acabado metálico cepillado. Se usa la foto."),
 "marco_puerto":     dict(d="Marcos de las ventanas/puertos del yugo", tipo="pintura",
                          foto="v6_p7_img0.jpeg:marco_ventana_foto", render="v6_p1_img0.jpeg:marco_ventana",
                          rough=0.35, fil="PLA gris azulado claro", conf="media",
                          nota="Foto: gris azulado claro. Render: gris oscuro (sombra en el render)."),
 "gris_azulado":     dict(d="Columnas de los gatos de elevación (pintadas en la máquina real)", tipo="pintura",
                          foto="v6_p7_img0.jpeg:rail_frontal_azulgris", render=None, rough=0.4,
                          fil="PLA gris azulado (RAL 7031/5014)", conf="media",
                          nota="En el render oficial son de acero; en la foto real, pintadas gris azulado."),
 "gris_azulado_claro": dict(d="Tapas/bombas de la parte superior (foto real)", tipo="pintura",
                          foto="v6_p7_img0.jpeg:bomba_superior_gris", render=None, rough=0.4,
                          fil="PLA gris azulado claro", conf="baja",
                          nota="Superficie curva con fuerte viraje azul en la foto: el tono real puede ser más gris."),
 "plastico_negro":   dict(d="Faldón inferior (lamas) y cadena portacables", tipo="plastico",
                          foto="v6_p7_img0.jpeg:cadena_cables_negra", render="v6_p1_img0.jpeg:faldon_negro",
                          rough=0.45, fil="PLA/PETG negro", conf="alta", pref="render",
                          nota="Se usa el render: las sombras de la foto real tienen un viraje azul de etalonaje."),
 "acero_cepillado":  dict(d="Patas, gatos (render), orejas de izado", tipo="metal", pbr="hierro",
                          foto=None, render="v6_p1_img0.jpeg:pata_acero", rough=0.33, fil="PLA seda plata", conf="media"),
 "inox":             dict(d="Cámara de vacío, carcasas RF (render), bridas, bombas, strippers", tipo="metal", pbr="inox",
                          foto="k2_p4_img0.jpeg:anillo_camara_inox", render="v6_p1_img0.jpeg:carcasa_rf_inox",
                          rough=0.18, fil="PLA seda plata", conf="alta", pref="render",
                          nota="La foto interior refleja el entorno azulado (variante 'foto'); color plano = render (neutro)."),
 "cobre":            dict(d="Dees y vástagos (cobre OFE), cátodos/chimeneas de fuente", tipo="metal", pbr="cobre",
                          foto="k2_p4_img0.jpeg:dee_cobre", render=None, rough=0.25, fil="PLA seda cobre", conf="alta"),
 "laton":            dict(d="Colector de agua de refrigeración (latón)", tipo="metal", pbr="laton",
                          foto="v6_p7_img0.jpeg:colector_laton", render=None, rough=0.3, fil="PLA seda oro/latón", conf="alta"),
 "acero_mecanizado": dict(d="Polos, sectores, tapón central (interior, acero 1010 mecanizado)", tipo="metal", pbr="hierro",
                          srgb=(150, 152, 156), rough=0.22, fil="PLA gris plata", conf="baja",
                          nota="No visible en las fichas (interior): acero mecanizado sin pintar, estimado."),
 "cobre_bobina":     dict(d="Bobina principal (cobre + resina, envuelta)", tipo="otro", srgb=(150, 85, 50), lin=(0.55, 0.28, 0.16),
                          metal=0.6, rough=0.45, fil="PLA cobre mate", conf="baja",
                          nota="No visible en las fichas: estimado (conductor de cobre impregnado)."),
 "gris_oscuro":      dict(d="Actuador de la fuente de iones", tipo="pintura", srgb=(55, 58, 62), lin=(0.16, 0.17, 0.18),
                          metal=0.3, rough=0.4, fil="PLA gris oscuro", conf="baja", nota="Estimado."),
 "niobio":           dict(d="Cuerpo de los blancos (Nb)", tipo="metal", pbr="niobio", srgb=(160, 162, 170), rough=0.3,
                          fil="PLA gris plata", conf="baja", nota="No visible en las fichas: niobio pulido, estimado."),
 "aluminio":         dict(d="Colimadores de los blancos", tipo="metal", pbr="aluminio", srgb=(200, 202, 206), rough=0.3,
                          fil="PLA seda plata", conf="baja", nota="Estimado."),
 "cuw":              dict(d="Chimeneas de la región central (CuW)", tipo="metal", srgb=(140, 110, 90), lin=(0.45, 0.33, 0.25),
                          metal=1.0, rough=0.35, fil="PLA bronce", conf="baja", nota="Estimado."),
 # servicios de la foto real (kiube/auxiliares.py)
 "manguera_azul":    dict(d="Mangueras de agua y neumática en espiral", tipo="plastico", foto="v6_p7_img0.jpeg:manguera_azul",
                          render=None, rough=0.5, fil="PETG/PLA azul", conf="media"),
 "valvula_roja":     dict(d="Manillas de las válvulas de bola del colector", tipo="pintura", foto="v6_p7_img0.jpeg:valvula_roja",
                          render=None, rough=0.4, fil="PLA rojo", conf="media",
                          nota="La corrección del viraje azul lo desplaza hacia naranja; el rojo real puede ser más puro (RAL 3000/3020)."),
 "transmisor_azul":  dict(d="Transmisores / caudalímetros de los ramales", tipo="plastico", foto="v6_p7_img0.jpeg:transmisor_azul",
                          render=None, rough=0.4, fil="PLA azul", conf="media"),
 "caja_electrica_gris": dict(d="Caja eléctrica (armario de chapa)", tipo="pintura", foto="v6_p7_img0.jpeg:caja_electrica_gris",
                          render=None, rough=0.4, fil="PLA gris claro / blanco", conf="media",
                          nota="Medido casi blanco; los armarios industriales suelen ser RAL 7035 (gris luminoso)."),
 "hormigon":         dict(d="Autoblindaje", tipo="otro", srgb=(170, 168, 160), lin=(0.62, 0.61, 0.58), rough=0.85,
                          fil="PLA gris piedra", conf="baja", nota="Las fichas no muestran el autoblindaje: estimado."),
}
EXTRA = {   # visibles en la foto pero no modelados (informativo, para pintar a mano o futuros LOD)
 "caja_azul_interior": "k2_p4_img0.jpeg:caja_azul_interior", "superficie_verdosa_interior": "k2_p4_img0.jpeg:superficie_verdosa",
}


def s2l(c):
    c = np.asarray(c, float) / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def l2s(c):
    c = np.clip(np.asarray(c, float), 0, 1)
    return np.round(255 * np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)).astype(int)


def lab(c):
    l = s2l(c)
    X = l @ np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]]).T
    X = X / np.array([0.95047, 1.0, 1.08883])
    f = np.where(X > 216 / 24389, np.cbrt(X), (24389 / 27 * X + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def ral(c):
    L = lab(c)
    k, v = min(RAL.items(), key=lambda kv: np.linalg.norm(lab(kv[1]) - L))
    return k, round(float(np.linalg.norm(lab(v) - L)), 1)


_img = {}
CMYK_ICC = "/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc"
SRGB_ICC = "/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/sRGB.icc"


def imagen(f):
    """sRGB; las imágenes CMYK de los PDF se convierten con el perfil FOGRA39 (impresión europea estándar)."""
    if f not in _img:
        im = Image.open(DOCS / f)
        if im.mode == "CMYK":
            from PIL import ImageCms
            im = ImageCms.profileToProfile(im, CMYK_ICC, SRGB_ICC, outputMode="RGB", renderingIntent=1)
        _img[f] = np.asarray(im.convert("RGB"), float)
    return _img[f]


def muestra_raw(f, reg):
    x0, y0, x1, y1 = R[f][reg]
    px = imagen(f)[y0:y1, x0:x1].reshape(-1, 3)
    if reg in CROMA:                                   # mitad más cromática
        ch = px.max(1) - px.min(1)
        px = px[ch >= np.median(ch)]
    lum = px @ [0.2126, 0.7152, 0.0722]
    lo, hi = np.percentile(lum, [25, 75])
    return np.median(px[(lum >= lo) & (lum <= hi)], axis=0)


def corregir(raw, f):
    """Corrección de 2 puntos en luz lineal: el blanco y el negro de referencia de la imagen pasan a neutros."""
    w_reg, k_reg = REF[f]
    lw = s2l(muestra_raw(f, w_reg)); lk = s2l(muestra_raw(f, k_reg)) if k_reg else np.zeros(3)
    Lw, Lk = lw.mean(), lk.mean()
    lin = Lk + (s2l(raw) - lk) * (Lw - Lk) / (lw - lk)
    return l2s(np.clip(lin, 0, 1)).astype(float)


def muestra(ref):
    f, reg = ref.split(":")
    raw = muestra_raw(f, reg)
    return raw, corregir(raw, f), f, reg


def albedo_pintura(bal, f):
    """lineal, con el blanco de referencia de esa imagen → 0.80"""
    w = s2l(corregir(muestra_raw(f, REF[f][0]), f)).mean()
    return np.clip(s2l(bal) * 0.80 / w, 0, 0.9)


def entrada(k, s):
    e = dict(descripcion=s["d"], tipo=s["tipo"], confianza=s["conf"], filamento=s["fil"])
    var = {}
    for v in ("foto", "render"):
        if s.get(v):
            raw, bal, f, reg = muestra(s[v])
            var[v] = dict(srgb=[int(round(x)) for x in bal], srgb_sin_balance=[int(round(x)) for x in raw],
                          hex="#%02X%02X%02X" % tuple(int(round(x)) for x in bal),
                          fuente=f"{FUENTE[f]} — región '{reg}' {list(R[f][reg])}")
            if not s.get("pbr") and not s.get("lin"):
                var[v]["base_lineal"] = [round(float(x), 3) for x in albedo_pintura(bal, f)]
    # color plano principal (impresión/visor): foto > render > estimado
    pref = s.get("pref", "foto")
    if pref in var:
        base = var[pref]
    elif "foto" in var:
        base = var["foto"]
    elif "render" in var:
        base = var["render"]
    else:
        base = dict(srgb=list(s["srgb"]), hex="#%02X%02X%02X" % tuple(s["srgb"]), fuente="estimado (no visible en las fichas)")
    e["srgb"] = base["srgb"]; e["hex"] = base["hex"]; e["fuente"] = base["fuente"]
    e["variantes"] = var
    r, de = ral(e["srgb"]); e["ral_aprox"] = f"RAL {r}"; e["ral_dE76"] = de
    # PBR
    if s.get("pbr"):
        e["pbr"] = dict(base_lineal=list(F0[s["pbr"]]), metalico=1.0, rugosidad=s["rough"],
                        fuente_base=f"F0 tabulado ({s['pbr']})")
    elif s.get("lin"):
        e["pbr"] = dict(base_lineal=list(s["lin"]), metalico=s.get("metal", 0.0), rugosidad=s["rough"], fuente_base="estimado")
    else:
        v = s.get("pref", "foto") if s.get(s.get("pref", "foto")) else ("foto" if s.get("foto") else "render")
        src = var[v]; f = s[v].split(":")[0]
        e["pbr"] = dict(base_lineal=[round(float(x), 3) for x in albedo_pintura(np.array(src["srgb"], float), f)],
                        metalico=0.0, rugosidad=s["rough"], fuente_base="muestra linealizada (blanco del yugo = 0.80)")
    if s.get("coat"):
        e["pbr"]["barniz"] = s["coat"]
    if s.get("nota"):
        e["nota"] = s["nota"]
    return e


pal = dict(
    version=1,
    fuentes={k: v for k, v in FUENTE.items()},
    metodo=__doc__.split("Uso:")[0].strip(),
    advertencias=["Las fichas técnicas no especifican colores ni códigos RAL: todo color es MEDIDO en imágenes "
                  "(fotos/render con iluminación, compresión JPEG y retoque de marketing) o ESTIMADO.",
                  "RAL más cercano: orientativo (ΔE76 contra tabla aproximada). Para pintar una maqueta, comparar con carta física.",
                  "Render oficial y foto real difieren (chaflanes, gatos, marcos, parte superior): variante por defecto = foto real."],
    materiales={k: entrada(k, s) for k, s in M.items()},
    extra_foto={k: (lambda r: dict(srgb=[int(round(x)) for x in r[1]], hex="#%02X%02X%02X" % tuple(int(round(x)) for x in r[1]),
                                   fuente=f"{FUENTE[r[2]]} — región '{r[3]}'"))(muestra(v)) for k, v in EXTRA.items()},
)
for k, v in pal["extra_foto"].items():
    r, de = ral(v["srgb"]); v["ral_aprox"] = f"RAL {r}"; v["ral_dE76"] = de
class _D(yaml.SafeDumper):
    def ignore_aliases(self, data):
        return True


yaml.dump(pal, Dumper=_D, stream=open(ROOT / "palette.yaml", "w"), allow_unicode=True, sort_keys=False, width=140)

# ---------------------------------------------------------------- figuras de control
(ROOT / "docs").mkdir(exist_ok=True)
for f, regs in R.items():
    im = Image.open(DOCS / f).convert("RGB"); d = ImageDraw.Draw(im)
    for k, b in regs.items():
        d.rectangle(b, outline=(255, 0, 0), width=2); d.text((b[0], max(0, b[1] - 11)), k, fill=(255, 0, 0))
    im.save(ROOT / "docs" / f"paleta_regiones_{f.split('.')[0]}.png")
items = list(pal["materiales"].items()) + list(pal["extra_foto"].items())
W, H = 1500, 46 * len(items) + 60
im = Image.new("RGB", (W, H), (250, 250, 250)); d = ImageDraw.Draw(im)
try:
    fn = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15)
except Exception:
    fn = ImageFont.load_default()
d.text((10, 15), "KIUBE: paleta medida en las fichas técnicas IBA (sRGB balanceado | PBR lineal→sRGB | RAL aprox.)", fill=(0, 0, 0), font=fn)
for i, (k, v) in enumerate(items):
    y = 50 + 46 * i
    d.rectangle((10, y, 110, y + 38), fill=tuple(v["srgb"]), outline=(0, 0, 0))
    if "pbr" in v:
        d.rectangle((115, y, 165, y + 38), fill=tuple(int(x) for x in l2s(v["pbr"]["base_lineal"])), outline=(0, 0, 0))
    if "ral_aprox" in v:
        d.rectangle((170, y, 200, y + 38), fill=RAL[v["ral_aprox"][4:]], outline=(0, 0, 0))
    conf = v.get("confianza", "foto")
    d.text((210, y + 2), f"{k}   {v['hex']}  {v['srgb']}   {v['ral_aprox']} (ΔE {v['ral_dE76']})   [{conf}]", fill=(0, 0, 0), font=fn)
    d.text((210, y + 21), v.get("fuente", "")[:150], fill=(90, 90, 90), font=fn)
im.save(ROOT / "docs" / "paleta_kiube.png")
print(json.dumps({k: (v["hex"], v["ral_aprox"], v["ral_dE76"]) for k, v in items}, indent=0, ensure_ascii=False))
