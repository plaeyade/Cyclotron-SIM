"""Regiones muestreadas (x0,y0,x1,y1 en píxeles) de las imágenes de las fichas técnicas IBA."""
R = {
 "v6_p1_img0.jpeg": {   # render oficial KIUBE (ficha VE24, p.1)
   "yugo_blanco_iluminado": (470, 450, 600, 520),
   "yugo_blanco_sombra":    (820, 450, 900, 520),
   "marco_ventana":         (279, 380, 291, 440),
   "faldon_negro":          (520, 585, 600, 600),
   "pata_acero":            (355, 640, 380, 700),
   "carcasa_rf_inox":       (622, 200, 650, 245),
   "banda_chaflan_render":  (480, 274, 700, 280),
 },
 "v6_p7_img0.jpeg": {   # foto máquina real (ficha VE24, p.7)
   "yugo_blanco_foto":      (380, 205, 450, 250),
   "banda_chaflan_foto":    (300, 163, 420, 170),
   "rail_lateral_azulgris": (127, 230, 143, 400),
   "rail_frontal_azulgris": (512, 300, 540, 400),
   "cadena_cables_negra":   (300, 148, 420, 158),
   "manguera_azul":         (296, 520, 312, 600),
   "colector_laton":        (335, 421, 395, 429),
   "valvula_roja":          (332, 447, 345, 455),
   "transmisor_azul":       (377, 560, 390, 585),
   "marco_ventana_foto":    (223, 303, 257, 317),
   "abertura_puerto_foto":  (193, 283, 213, 297),
   "caja_electrica_gris":   (200, 490, 240, 560),
   "bomba_superior_gris":   (345, 112, 385, 128),
 },
 "k2_p4_img0.jpeg": {   # foto interior (ficha KEY22, p.4)
   "dee_cobre":             (500, 1092, 700, 1125),
   "yugo_blanco_interior":  (150, 450, 330, 500),
   "agujero_camara":        (190, 872, 215, 892),
   "superficie_verdosa":    (20, 1060, 140, 1100),
   "anillo_camara_inox":    (300, 818, 500, 836),
   "caja_azul_interior":    (70, 600, 220, 700),
   "asa_roja":              (389, 565, 399, 630),
 },
}

# Referencias neutras por imagen para la corrección de 2 puntos (blanco, negro) en luz lineal
REF = {"v6_p1_img0.jpeg": ("yugo_blanco_iluminado", None),
       "v6_p7_img0.jpeg": ("yugo_blanco_foto", "abertura_puerto_foto"),
       "k2_p4_img0.jpeg": ("yugo_blanco_interior", "agujero_camara")}
# Regiones metálicas/saturadas: se toma la mitad más cromática (evita reflejos blancos del entorno)
CROMA = {"asa_roja", "dee_cobre", "colector_laton", "manguera_azul", "valvula_roja", "transmisor_azul"}
