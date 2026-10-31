# Modelos de elementos finitos (equivalentes a Opera-3D) y región central

Opera-3D (Simulia) es comercial y no corre aquí. Cada modelo se entrega dos veces:

- como **paquete para Opera-3D**: geometría STEP, materiales, conductores o potenciales, contorno, malla y exportación;
- como **gemelo ejecutado** con NGSolve (FE de código abierto, misma física).

El gemelo trae resultados de referencia y un comparador, para verificar el modelo de Opera en cuanto lo corras.

| Modelo | Gemelo ejecutado | Paquete Opera-3D | Comparador |
|---|---|---|---|
| Imán del KIUBE (magnetostático no lineal, sector de 90° periódico, z ≥ 0) | `fem/magnet_ngsolve.py` (Nédélec, Newton con curva B-H) | `out/opera3d/OPERA_SETUP.md` | `fem/opera_compare.py tabla_opera.csv` |
| Región central C18/9 (electrostático 3D, malla graduada) | `fem/central_region.py --preset c18_9` (H1 de orden 2) | `out/opera3d/region_central_c18_9/OPERA_REGION_CENTRAL.md` | `fem/central_region.py --campoE tabla_opera.txt` |

## Imán

```bash
python build.py --lod fem_mag                                   # dominios STEP (sector 90° y mitad de 360°)
python fem/magnet_ngsolve.py --diag                             # NI_FE_calibrado (116 kA·v), acero del yugo, reparto de Ampère
python fem/magnet_ngsolve.py --NI 105000 116000 --curva bh_steel_1010 --diag   # comparar con la 1010 genérica
python fem/fem_mag_variante.py --tag prueba --set pole.valley_floor_z=100       # estudios de sensibilidad
python validation/calibrar_iman.py                                # tabla y figura de calibración frente a KEY22/S09
python fem/opera_export.py                                        # paquete Opera-3D/TOSCA
```

### Calibración final del imán (v0.21): 19 kW, patentes y KEY22

Decisión del usuario: la potencia de referencia es la de IPAC26, **19 kW de bobinas**. Con 2 t de Cu, una vuelta media de 4.17 m y ρ(40 °C) eso da NI = 115 kA·v. Las patentes de IBA de 2016 (`PAT16`) fijan dos cosas que v0.18 no cumplía:
- **gap de valle Gv = 200–250 mm** (US9907153, reivindicación 19);
- **corrector de gradiente en cuña de 90°** con abertura L10 = 3–3.5·Gh y profundidad H10 = 6–10 % de Lh (US10278277). Ahora es `pole.gc_forma: cuna`, con 78 × 39 mm en la mitad de la colina del lado del gap.

| Caso a J_s = 2.12 T, tapa 181 | NI | ⟨B⟩ 0.38–0.47 m | Flutter | Colina / valle (0.45 m) | ⟨B⟩(0)/⟨B⟩(0.45) |
|---|---|---|---|---|---|
| v0.18: valle 95 (Gv 190), ranura 48 × 26 | 116 k | 1.3240 | 0.381 | 2.20 / 0.446 | 1.026 |
| Valle 110 (Gv 220) + cuña | 115 k / 125 k | 1.3102 / 1.3365 | 0.397 / 0.384 | 2.20 / 0.410 | 1.009 |
| Valle 100 (Gv 200) + cuña | 115 k | 1.3250 | 0.390 | 2.22 / 0.425 | 1.011 |
| **v0.21: valle 100 + cuña + tapón −15** | **115 k (= 19 kW)** | **1.3253** | **0.389** | **2.22 / 0.426** | **1.030** |
| Malla de 6 mm (mismo caso) | 115 k | 1.3246 (−0.05 %) | 0.390 | 2.22 / 0.423 | — |
| **KEY22 / S09** | — | **1.3253** | **0.394** | **2.15 / 0.45** | **1.033** |

El modelo **cierra sin tocar el acero**: con las cotas de las patentes, el campo de IBA sale exactamente con la corriente que permiten 19 kW. La cuña mejora el flutter cerca de la extracción (0.45–0.5 m; `validation/calibracion_iman.png`). Gv = 200 está en el borde inferior de la reivindicación; con Gv = 220 harían falta 121 kA·v (21 kW). La discrepancia con «< 11 kW» de TUD03T sigue registrada en `docs/FUENTES.md`.

### Calibración del circuito magnético (v0.18)

Hasta v0.17 el FE pedía 170 kA·v para el campo de IBA, y la potencia de las bobinas solo da unos 105 kA·v. El diagnóstico `--diag` reparte ∮H·dl = NI/2 a lo largo del camino del flujo: gap → polo → sector C → tapa A → yugo lateral B → plano medio. El cierre es mejor que el 1 %.

| Tramo (mitad superior) | 1010 genérica, 105 kA·v | 1010 genérica, 170 kA·v | Calibrado, 116 kA·v |
|---|---|---|---|
| Gap de colina (aire) | 19.8 kA | 21.4 kA | 21.7 kA |
| Polo + sector C (hierro) | 17.7 kA (2.21 T) | 30.8 kA (2.35 T) | 20.1 kA (2.39 T) |
| Tapa A, radial | 14.8 kA (2.12 T) | 31.9 kA (2.22 T) | 15.9 kA (2.24 T) |
| Yugo lateral B | 0.4 kA | 1.2 kA | 0.7 kA |
| Fracción en el hierro | 63 % | 75 % | 63 % |

Las puntas de polo y la tapa trabajan por encima de 2.1 T, donde manda la polarización de saturación J_s. La curva 1010 genérica satura en J_s ≈ 2.0 T, y el acero de bajo carbono de los ciclotrones llega a 2.10–2.15 T (hierro puro: 2.16 T).

**Sensibilidades a 105 kA·v** (`out/fem_mag/calib/`, métricas en la banda 0.38–0.47 m como KEY22):

| Caso | ⟨B⟩ [T] | Flutter |
|---|---|---|
| 1010 genérica (J_s 2.00) | 1.187 | 0.415 |
| J_s 2.12 | 1.242 | 0.423 |
| J_s 2.12 + tapa 181 | 1.271 | 0.409 |
| J_s 2.12 + valle 90 | 1.269 | 0.407 |
| J_s 2.12 + tapa 181 + valle 95 | 1.293 | 0.396 |
| **KEY22** | **1.325** | **0.394** |

Por sí solo J_s no basta: haría falta unos 2.3 T, por encima del hierro puro. La calibración combina tres palancas físicas:

- **Acero:** `bh_acero_bajoC_kiube` en `materials.yaml`, con J_s = 2.12 T y la misma forma que la 1010.
- **Tapa A:** 181 mm, el extremo superior de 166 ± 15 de T04F1.
- **Valle:** a 95 mm, fijado por el flutter. Con la curva genérica salía 110 porque compensaba la saturación de más.

Con el acero corregido, el tapón pasa de −17 a −16 mm. El NI se interpola para ⟨B⟩ = 1.3253 T y sale **116 kA·v**.

| Resultado a 116 kA·v | FE v0.18 | Referencia |
|---|---|---|
| ⟨B⟩ en 0.38–0.47 m | 1.3240 T | 1.3253 T (KEY22) |
| Flutter en 0.38–0.47 m | 0.381 | 0.394 (KEY22) |
| B colina / B valle en 0.45 m | 2.20 / 0.446 T | 2.15 / 0.45 T (S09) |
| ⟨B⟩(0)/⟨B⟩(0.45) | 1.026 | 1.033 (KEY22) |

Los 116 kA·v quedan un 10 % por encima de los 105 de la potencia, que es una estimación con un rango de 100 a 115. La diferencia se reduce con una bobina algo más alta (190 mm es EST; con 205 mm la misma potencia da 109 kA·v) o con un acero de J_s = 2.15 T (que baja la necesidad a unos 111 kA·v). Juntos la cierran.

**La calibración no es única:** J_s, el espesor de la tapa y la altura de la bobina se compensan entre sí. Una curva B-H medida del lote fijaría J_s y dejaría solo la geometría.

**Chequeo CY11** (Cyclone 11 de IBA, TOSCA): su campo sube 160 G por 10 A a 186 A, es decir d ln⟨B⟩/d ln NI = 0.21. El FE calibrado del KIUBE da 0.24 entre 105 y 115 kA·v: un nivel de saturación del hierro parecido al de una máquina de IBA. En el C11, rellenar el valle 50 mm cambió el flutter menos de un 5 %; aquí, subir el valle 15 mm lo cambia un −3 %. IBA usa la curva B-H del proveedor del hierro: es el dato que cerraría la no unicidad.

**Chequeo BATAN13** (Opera-3D/TOSCA, 13 MeV): B_valle/(μ0·NI/g_valle) vale 0.84 en esa máquina. En el KIUBE calibrado da 0.58, porque su valle profundo (190 mm) recibe el campo sobre todo de las colinas por el borde y no de la bobina directamente. La comparación sirve como cota y no como objetivo.

La Ec. 31 de CAS 2017 (`validation/caida_hierro_ec31.json`) da para el C18/9 una caída en el hierro del 60 %. El FE calibrado del KIUBE da un 63 %, coherente.

**Opera-3D/TOSCA:** es comercial y no corre aquí. Su paquete queda listo en `out/opera3d/` con la curva calibrada y NI = 116 kA·v, y NGSolve lo sustituye mientras tanto.

## Región central (Cyclone 18/9 y KIUBE)

**v0.22: preset `kiube` con el mapa Opera de IBA** (TUD03T, diap. 5): dees de 50° cerca del centro, puntas a r = 8 mm y fuente a r = 12.5 mm. La fuente queda 39.4° antes del borde de entrada de su dee (IBA: 40°). Para mallas de más de 250 k tetraedros el campo se resuelve con CG + h1amg, porque bddc no cabe en 7 GB.

| Fuente a | Descentrado tras 6 vueltas | Fases que sobreviven |
|---|---|---|
| 10.0 mm | 11.4 mm | 7 |
| 11.25 mm | 8.5 mm | 7 |
| **12.5 mm (IBA)** | **6.2 mm** | 5 |
| 13.75 mm | 4.8 mm | 7 |
| 15.0 mm | 4.4 mm | 7 |

**v0.25: rendija de la chimenea.** Con el puller final, la pérdida vertical ocurre casi toda en la primera vuelta (r = 20–30 mm) y depende de la altura de salida del ion: sobrevive el 49 % con |z₀| < 0.8 mm, el 18 % con 0.8–1.6 mm y el 9 % con más de 1.6 mm. Por eso se barrió la altura de la rendija (`validation/region_central_rendija.png` y `.json`):

| Rendija | 3 | 4 | 5 | **6** | 8 mm |
|---|---|---|---|---|---|
| Transmisión (± 2.5 %) | 54.8 | 43.8 | 37.2 | **34.8** | 25.0 % |
| Corriente útil ∝ altura × transmisión | 1.64 | 1.75 | 1.86 | **2.09** | 2.00 |

Una rendija baja transmite mejor, pero la fuente PIG emite en proporción al área de la rendija. **La corriente que llega a la extracción es máxima con 6–8 mm**, así que la rendija se queda en 6 mm (CALC). La transmisión del ≈ 35 % es el precio de esa corriente, no un defecto del puller.

**v0.24: forma del puller.** Barrido con el gap de 2.8 mm, malla fina, 6 vueltas y 400 partículas (`fem/central_region_barrido_puller.py`; figura `validation/region_central_puller.png`, datos `.json`).

| Abertura vertical | 3 | 4 | 5 (v0.23) | 6 | 7 | 8 | **10** | 12 mm |
|---|---|---|---|---|---|---|---|---|
| Transmisión (± 2 %) | 14.5 | 16.8 | 19.5 | 19.8 | 26.2 | 26.0 | **34.8** | 31.0 % |
| Choca con metal | 72 | 64 | 53 | 40 | 30 | 17 | **3** | 1 % |
| Se pierde en vertical | 14 | 20 | 28 | 40 | 44 | 57 | **62** | 68 % |

Con la abertura de 5 mm, el 37 % del haz moría a r < 15 mm, en las mordazas del puller. Con 10 mm el puller casi no intercepta (3 %) y el límite pasa a ser vertical. Con la abertura fijada en 10 mm:
- **ancho:** 5 → 29.5 %, **7 → 34.8 %**, 10 → 33.5 %;
- **profundidad:** 2 → 33.0 %, **3 → 34.8 %**, 5 → 29.8 % (con 5 el centrado empeora a 13 mm).

**Abertura de los dees en el centro.** Con ±14 mm (la de la caja del dee en el CAD), el campo acelerador no alcanza la chimenea: la aceptancia cae a 20–30°, el 82–90 % del haz choca en la primera vuelta y la transmisión baja al 8–14 %. Por eso las puntas de los dees del CAD pasan a z = 8–12 (`dee.punta_z0` = 8, v0.24), con un escalón de unión a la caja (que sigue a ±14). La RF se recalibra: la caja pasa a 78.7 mm y da 40.652 MHz.

Puller final (CALC): abertura de 10 mm, ancho de 7 mm, profundidad de 3 mm y gap de 2.8 mm. La transmisión pasa del 19.5 al ≈ 35 %, frente a un techo de ≈ 87 % por aceptancia de fase. Lo que queda es pérdida vertical (62 %), que depende de la dispersión vertical inicial en la rendija (EST, ±2.4 mm) y del enfoque de las primeras vueltas.

**v0.23: gap fuente–puller.** Barrido con malla fina, 6 vueltas y 400 partículas, con la fuente en la posición de IBA (r = 12.5 mm, 40° antes del borde). Resultado en `validation/region_central_gap.png` y `.json`:

| Gap | Transmisión | Pérdidas en metal | Pérdidas verticales | Campo en el gap |
|---|---|---|---|---|
| 2.0 mm (campo del C18/9 escalado) | 14.8 ± 1.8 % | 62 % | 23 % | 20.0 MV/m (2.4 Kp) |
| 2.3 mm | 18.2 ± 1.9 % | 57 % | 25 % | 17.4 MV/m |
| 2.6 mm (v0.22) | 21.5 ± 2.1 % | 55 % | 24 % | 15.4 MV/m |
| **2.8 mm (v0.23)** | **19.5 ± 2.0 %** | **53 %** | **28 %** | **14.3 MV/m (1.7 Kp)** |
| 3.2 mm | 19.2 ± 2.0 % | 48 % | 33 % | 12.5 MV/m |

Hay un óptimo (ajuste parabólico con los 5 puntos: **2.79 mm**, 20.6 %). Con gaps cortos el haz choca contra el puller; con gaps largos el campo de extracción enfoca menos en vertical y crecen las pérdidas en z. **Corrección (v0.24):** las partículas del haz se lanzan con fases en la ventana aceptada ±10° (≈ 120°), no en 360°, así que el techo por aceptancia de fase es ≈ 105/120 ≈ 87 % y no el 29 % que se dio en v0.23. Una transmisión del ~20 % significa que se pierde mucho haz, y sobre todo en el paso fuente–puller: 37 % del haz termina a r < 15 mm (v0.24). Con 2.8 mm, la regla del modelo (puller pegado al borde de entrada) coloca la fuente a 40.3° del borde, igual que el mapa de IBA, sin ajuste extra. El gap del C18/9 (1.5 mm a 32 kV) escalado al campo de 40 kV (2.0 mm) transmite un 7 % menos.

El caso nominal de v0.22 (12.5 mm, malla fina) acepta 105° de fase RF y transmite el 18 % del haz en 6 vueltas. Con 15 mm, el 17 %. El modelo centra mejor con la fuente a 14–15 mm. Como la chimenea, la rendija, el gap fuente–puller (2.6 mm) y el puller son EST, la geometría conserva la posición medida de IBA; la diferencia apunta a esas cotas. Por ejemplo, el C18/9 usa un gap de 1.5 mm (JUAS22). El cuello que une los dees no está en este modelo electrostático: está al mismo potencial que los dees y cerca del eje.


```bash
python fem/central_region.py --preset c18_9 --optimizar 0.6 0.8 1.0 1.2 1.4   # optimiza la posición de las fuentes
python fem/central_region.py --preset c18_9 --set gap_sp=3.0 V_dee=30000      # cualquier cota es un parámetro
python fem/region_central_escena.py out/central_region/c18_9                  # escena 3D para Blender
blender -b -P blender/render_region_central.py -- --escena out/central_region/c18_9/escena/escena.json \
        --salida rc.png --res 3840x2160 --muestras 256 --gpu auto
```

**Física**

- **Campo eléctrico:** Laplace 3D con dees y pullers a 1 V, y a tierra las colinas, el tapón, los dummy dees y las chimeneas. Simetría de plano medio (E_z = 0 en z = 0).
- **Malla graduada:**
  - 0.2 mm en cada gap fuente–puller y su rendija;
  - 1 mm en las narices de los dees;
  - 5 mm hasta r = 150 mm;
  - 18 mm en el resto de la estructura de los dees, hasta r = 480 mm.
- **Campo magnético:** B0 isócrono por especie (H⁻ con h = 2, D⁻ con h = 4), con joroba central de 250 G y colinas de 1.9 / 0.35 T suavizadas. Br y Bφ salen de la expansión de Maxwell en z.
- **Tracking:** Lorentz relativista (Boris) con E(x, t) = E₁(x) · V_dee · cos(ωt + φ₀), evaluado directamente en el elemento.

**Datos públicos del C18/9** (lista de ciclotrones de Cyclotrons 2004, CM06, y 2001, CU-08):

- 42 MHz, con h = 2 para p y h = 4 para d;
- 32 kV, dees de 30°;
- ⟨B⟩ 1.35 T, colina 1.9 T y valle 0.35 T;
- gap de colina 30 mm, radio de extracción 0.48 m, 2 PIG internas.

La geometría fina (chimenea, rendija, puller, narices, dummy dees) **no es pública**: es diseño típico, parametrizado.

**Mejoras tomadas de Zaremba & Kleeven, CAS 2017 (arXiv:1804.08961, §3.3)**

- **Ajuste longitudinal (Fig. 37):** los puntos V = 0 deberían caer sobre los ejes (cada 180°/h). Se mide la desviación media y se reporta qué fracción de la ventana aceptada cae en la guía −40…−10°.
  - Un barrido del azimut de las fuentes fijó +4°. La desviación pasó de 12.9° a 7.6° (H⁻) y de 8.9° a 8.4° (D⁻), con fase de diseño D⁻ en −35°.
  - El detalle está en `validation/region_central_estudio_parametrico.json`.
- **Recorte de la chimenea D⁻ para el paso del H⁻ (Fig. 37),** parámetro `corte_ext`.
  - Sin recorte, el 13 % del haz H⁻ muere en la chimenea D⁻ (20 de 150 partículas). Con 1.6 mm baja a 15, sin cambio significativo de la transmisión total en esta geometría estimada. Queda en 0 por defecto.
- **Fuerza vertical eléctrica a lo largo de la órbita (Fig. 40):** `fuerza_vertical.png`. Cada gap enfoca a la entrada y desenfoca a la salida, y el efecto cae con las vueltas.
- **Chimeneas de Cu-W y cátodos de Ta (§3.3.1),** en el render.

**Figuras** (`out/central_region/c18_9/`):

- `orbitas_planta.png`: equivalente a CAS 2017 Fig. 37. Puntos rojos donde V = 0 y verdes donde V es máximo.
- `fuente_puller_zoom.png`, `aceptancia_energia.png`, `vertical.png`, `fuerza_vertical.png`.
- `region_central_c18_9_3d.png`: render de Blender.

### v0.29: nariz del dee, chimenea y diagnóstico vertical

- Preset `kiube` con la geometría del CAD:
  - `puente_w` 3.5 y `puente_off` 1.25 (cuello entre puntas);
  - `r_caja` 60, `ap_caja` 14 y `dee_angle_caja` 39 (caja del dee y escalón);
  - `plug_recess` 15 (tapón rebajado) y `H_dom` 32.
- El límite vertical de pérdida depende del radio: ±(ap_half − 0.2) dentro de las puntas y ±(min(ap_caja, hill_gap_half) − 0.2) fuera.
- `slit_rot` gira la rendija y el puller respecto de la tangente (+ = hacia fuera).
- Barridos: `python fem/central_region_barrido_puller.py <tag> clave=valor ...`. Figura: `python validation/figura_region_central_v29.py`.
- `python fem/central_region_diag_z.py` muestra dónde y cuándo se pierden las partículas por |z| (`out/central_region/diag_z/`).

### v0.30: menisco del plasma (fuente PIG)

- `python fem/menisco.py [--jp 1440] [--cos 0.97 0.86 0.7 0.5]` calcula el menisco autoconsistente por nivel de voltaje y escribe `out/central_region/menisco/menisco_jp<jp>.json`.
  - La primera vez resuelve el campo global y guarda la rejilla del potencial en `rejilla_phi_global.npz`; `--recalcular` la rehace.
  - El dominio local es la chimenea, el canal de la rendija y las mordazas del puller; en la caja se impone el potencial global.
  - El plasma se representa por penalización de volumen (φ = 0), así que la malla es fija.
  - La emisión es min(j_p, j_CL local), con carga espacial depositada por macropartículas en un campo B de 1.35 T.
- `python fem/central_region_menisco.py [--jp 1440] [--plano]` hace el seguimiento global con el menisco: un campo por nivel y las partículas repartidas por cos φ. `--plano` es el control con el plasma plano.
- En el preset: `menisco_s_cav` (fondo de la cavidad de la rendija, 1.8 mm) y `_menisco` = [c0, c1, c2, c3, a_u], con s(u, z) = c0 + c1 z² + c2 z⁴ + c3 z⁶ + a_u (2u/w)². Además, el canal de la rendija se malla a `h_canal` = 0.08 mm.
- **Limitaciones:** plasma frío; electrones extraídos ignorados; carga espacial solo en el dominio local; menisco cuasiestático por nivel de voltaje.

### v0.31: forma de las mordazas del puller

- `mordazas_puller(pr, z_top)` devuelve las 2 mordazas en el marco local (x = 0 en la cara frontal). La usan `build()` y `fem/menisco.py`.
  - `puller_h2` y `puller_d1`: abertura más estrecha o más ancha a partir de x = d1.
  - `puller_chaflan` y `puller_chaflan_ang`: bisel del borde interior frontal, con catetos c en x y c·tan(ang) en z.
- Preset `kiube`: bisel de 5 mm a 45°.
- `python fem/central_region_campo_max.py <tag> [clave=valor]` calcula |E| máximo en una rejilla de 0.05 mm entre la chimenea y el puller a 40 kV.
- Figura: `python validation/figura_region_central_v31.py`.

### v0.33: órbitas, isocronía y aceleración

- `python fem/orbitas.py` calcula las órbitas de equilibrio de 0.5 a 18.5 MeV en el mapa FE y escribe `out/orbitas/equilibrio.json`. `orbitas.mapa_key22()` construye el mapa sintético con el ⟨B⟩ y el flutter medidos (validación).
- `python fem/aceleracion.py --npz … --escala k --fases a b --n N` sigue el paquete en 3D hasta 18 MeV.
- `python fem/shim_iso.py respuesta i | caso TAG '[Δα…]' [pole.x=v] | ajustar BASE '[Δα…]' [fijas]` hace el ajuste de las colinas. Cada caso es un FE del sector de 90° (≈10 min) más sus órbitas; los resultados van a `out/fem_mag_variantes/<TAG>/`.
- La perilla de 60 mm se descarta porque deja una arista fina junto al tapón. Algunas combinaciones del borde dan elementos degenerados; el script reintenta con otra malla.

### v0.34: ajuste fino con jacobiano de órbitas

- `python fem/shim_iso2.py jac i` calcula el caso con la perilla i perturbada (1–10: Δα +0.7° por lado en `R2`; 11: L +4 mm) sobre la base `shim_it3_L520`. Cada caso es un FE de ≈10 min más sus órbitas.
- `resolver(b, cols, base_da, base_L, reg, peso_fase, regL, regc, cotas=cotas_dee())` resuelve el paso:
  - filas de error punto a punto (1 = 100 ppm) y de fase acumulada del ion central (1 = 20°; ΔE = 4·40 kV·sin 39° por vuelta);
  - incógnitas δ, nivel c y fase inicial;
  - Tikhonov sobre δα, δL y c, más segunda diferencia del perfil;
  - con `cotas`, `scipy.optimize.lsq_linear` con Δα ≤ máximo por holgura de 5 mm a la caja del dee y 515 ≤ L ≤ 526.
- Iteraciones en `out/fem_mag_variantes/shim_it4…it8`. it7 (sin cotas) invadía los dees; **it8** es la adoptada.
- El mapa oficial (`out/fem_mag/ngsolve`, dominios regenerados) reproduce it8 con σ = 0.07 % y máximo 0.13 %; `out/orbitas/equilibrio_v34.json`.

### v0.35: fuente → blanco

- `python fem/fuente_a_blanco.py rc --n 800 --escala 1.0103 --rc rc_k1.0103.npz` corre la región central con el menisco en el mapa FE oficial (4 cálculos FE, uno por nivel del menisco; ≈30 min). Guarda `out/fuente_blanco/rc_*.npz` con X, V, T, t, fase de emisión y motivo.
- `... acel --escala k --suf S --rc ...` sigue el haz hasta la lámina y escribe `acel_p0S.json` con `estados_lamina`. `... ext --escala k --suf S` hace el barrido de azimut de la lámina, el ajuste y la mancha en el blanco, y escribe `ext_p0S.json`.
- `central_region.py`: con `B_source="fem_map"`, las claves `B_npz` y `B_escala` fijan el mapa y el nivel (antes se normalizaba a B0 en 20–80 mm).
- `aceleracion.acelerar(..., estado=dict(X, P, T, phi0), lamina=dict(phi_deg, r_in, r_out), r_max=0.60)`.
- `magnet_ngsolve.py` escribe además `bz_midplane_ext_*.npz`: 0.40–0.95 m, puntos fuera de la malla → NaN, y el material.
