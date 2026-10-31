# Registro de fuentes documentales

> **Uso personal.** Todos los documentos son propiedad del usuario, que cuenta con los derechos de uso. El modelo es de uso personal y no público, por lo que no se restringe información confidencial en la construcción ni en las adecuaciones.

## Revisados en v0.27

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `VUB24` | Covens (VUB-UZ Brussel), *Target shielding in a PET-production facility*, FANC 16-10-2024 (19 diap.; `ed9fd30a…`) | **KIUBE real** (puesto en marcha en 2017): blancos cónicos Nirta 2XL (puertos 3 y 7) y 3XL (4 y 8); **puertas blindadas de HDPE de 18 cm** (ficha IBA); foto de las puertas en la cara del yugo; carga de 50 000 µAh/año en T4/T8; **activación del hormigón medida** (Eu-152, Co-60 a 0–100 cm) y calculada por IBA (÷2.3 con puertas) | **Incorporado:** puertas de HDPE modeladas (`target.puerta_hdpe_mm`) y datos de activación en `validation/activacion_vub_kiube.json` para validar el blindaje |
| `CON13` | Devillet et al. (IBA), *Performance of IBA new conical shaped niobium [18O] water targets*, Cyclotrons2013 WEPSH043 (`d5ca725c…`) | Insertos cónicos de 2.4/3.4/5/7 ml; **ventana Havar de 35 µm en operación** (no 50); hasta 145 µA; rotura > 50 bar | **Incorporado:** `tgt_window_um` = 35 |
| — | VE24 (ya registrado), Figs. 3 y 4, releídas | Fig. 3 (ejes en m): **los blancos están fuera de la cara del yugo** (±0.88–0.98 m); Fig. 4: cámara de vacío extendida por el plano medio hasta el yugo | **Incorporado:** `extraction.blanco_en_cara` (los blancos pasan de r = 590 a la cara del yugo, +10 mm), con tubo de haz por la escotadura |
| — | Tochon-Danguy et al., *Implementation of a solid target production facility* (AIP 2012); Gelbart & Johnson, *Solid target system with in-situ dissolution* (Instruments 2019); Dellepiane et al., *Novel solid target… Bern* (EPJ 2023) | Estaciones de blanco sólido (Nirta STS y otras), haz de 13–16 mm, discos de Ø24 mm | **Referencia** (el proyecto no incluye blanco sólido) |
| — | Sumitomo, *SC230* (AFAD 2024); Smirnov & Vorozhtsov, *Modern compact accelerators of cyclotron type* (PPN 2016) | Ciclotrones superconductores de protonterapia y compactos | **Desechado** (otras máquinas) |

**Benchmark Geant4 frente al SRG, repetido en v0.27** con la fuente en la cara del yugo (r = 969 mm en lugar de 597), sin puertas, como las Tablas 8–9:
- dentro del búnker el modelo pasa de 0.23–0.54 a **0.51–1.22 del SRG**: buena parte de la diferencia anterior era la posición del blanco;
- tras 200–230 cm queda a 0.06–0.08 del SRG (antes 0.007–0.02); el TVL del modelo es de 28.6 cm (ajuste ponderado) frente a 33.9 cm del SRG;
- **con las puertas de HDPE de 18 cm** la fluencia de neutrones a 0–10 cm del muro frontal baja ÷1.85–2.0 (IBA: activación ÷2.3, VUB24) y el perfil en profundidad (1 / 0.45 / 0.13 / 0.035 por tramos de 15 cm) cae como las varillas medidas en Bruselas (1 / 0.41–0.57 / 0.09–0.19 / 0.07) → `validation/benchmark_puertas.{json,png}`;
- detalle en `mc_benchmark/README_BENCHMARK_SRG.md`; resultado anterior en `validation/benchmark_srg_v026_r597.json`.

## Menisco de la fuente (v0.30)

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `KUO92` | Kuo (TRIUMF) y Hendry (Cyclotron Inc.), *Recent development of internal H⁻/D⁻ source for compact cyclotrons*, Cyclotrons 1992, pp. 352–356 ([KEK](https://epaper.kek.jp/c92/papers/viii-13.pdf)) | Fuente H⁻ interna con rendija de 0.6 × 4 mm: **2.8 mA equivalentes DC, 144 mA/cm²**; gap de expansión del plasma de 0.6–1 mm; ancho de fase de unos 64° para 500 µA a r = 15 cm (32 kV) | **Incorporado:** j_p de referencia del menisco (`fem/menisco.py`); sensibilidad ×0.5 / ×2 |

## Techo del búnker (v0.28, preliminar)

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `RGSR` | *Reglamento General de Seguridad Radiológica* (México, DOF 22-11-1988), [texto en tabasco.gob.mx](https://tabasco.gob.mx/leyes/federales/reglamento/pdf/1526) | POE 50 mSv/año (Art. 20); público 1/10 (Art. 37: 5 mSv/año) y 1 mSv/año de promedio en exposición prolongada (Art. 38); zona de radiación > 0.05 mSv/h o 1 mSv en 5 días (Art. 113) | **Criterio de diseño:** 1 mSv/año en el cuarto de A/A y < 50 µSv/h durante la irradiación |
| `NCRP151` | NCRP Report 151 (2005), tabla 4.1 de factores de ocupación (vía [SAHPRA, *Shielding*](https://www.sahpra.org.za/wp-content/uploads/2019/09/Shielding.pdf)) | T = 1/40 en áreas de ocupación mínima (áticos, escaleras); 1/20 en almacenes | Cuarto de A/A sin puesto de trabajo: **T = 1/40** (1/20 como variante) |

**Cálculo** (`mc_benchmark/bunker_techo`, `validation/techo_bunker.{json,png}`):
- búnker del proyecto 5.9 × 4.7 × 2.5 m, ciclotrón centrado, blancos T1/T5 con puertas de HDPE, hormigón 2.35 g/cm³;
- TVL del techo de 24.6 cm;
- espesor para 1 mSv/año con 50–150 k µAh/año: 1.27–1.46 m (modelo) y **1.55–1.74 m escalado al SRG** (×14, el sesgo del modelo tras 2 m en el benchmark);
- < 50 µSv/h a 2 × 100 µA: 1.68 m (escalado).

## Patentes de IBA de 2016 (v0.20b)

Son las 4 solicitudes que cita TUD03T, todas con prioridad del 13-05-2016 y publicadas como EP3244707–710. Se leyeron sus equivalentes en EE. UU. **Ninguna da cotas absolutas del KIUBE, ni potencia ni NI**: solo rangos preferidos y reivindicaciones pensados alrededor de una máquina de 18 MeV como el KIUBE. La comprobación está en `validation/patentes_iba.py` → `.json`.

| Prioridad | Publicación EE. UU. | Título | Qué aporta |
|---|---|---|---|
| EP16169489 | US10064264 | Pole insert for cyclotron | Inserto de 15–150 mm de ancho (≤ 15 % de Ah), altura 40–60 % de la colina, largo 400–800 mm, mismo acero que la colina, atornillado y mecanizado iterativamente en el mapeo. El modelo (52 × 35 mm) cumple |
| EP16169490 | US9961757 (`PAT`, ya registrada) | Peripheral hill sector design | Rh/Lh ≤ 65 %; borde de la colina concéntrico con la órbita, a ≈ 0.6·Gh de ella. El modelo tiene Rh/Lh = 0.64 y el centro del arco en el centro de curvatura de la órbita (coincide exactamente). Separación borde–órbita: 41 mm en el modelo de borde duro frente a 14 mm de la regla, lo que sugiere que la órbita real de 18 MeV en la colina está a ≈ 498 mm |
| EP16169494 | US10278277 | Gradient corrector for cyclotron | Corrector rebajado en cuña de 90°: profundidad 6–10 % de Lh, ancho 3–3.5·Gh, 45–55 % de la altura; 2 por colina (8 salidas); reduce la potencia de la bobina (cualitativo). El modelo (26 × 48 mm) queda por debajo |
| EP16169497 | US9907153 | Compact cyclotron («geometry») | Gh 22–26 mm; **Gv 200–250 mm (reiv. 19)**; αv 42–45°; agujeros de bombeo 2Ra 240–300 mm, 2Ra/La ≥ 70 % (reiv. 16); (Gv·Tv)/Lv² < 10 %; (Gh·Gv)/Lv² < 2 %; Lh 520–550 mm; ≈ ⅓ más ligero que la generación anterior |

**Aplicado en v0.21** (con 19 kW de IPAC26 como referencia, decisión del usuario): valle a 100 mm (Gv = 200) y corrector de gradiente en cuña de 78 × 39 mm. El FE da ⟨B⟩ = 1.3253 T y flutter 0.389 a 115 kA·v = 19 kW, sin cambiar el acero. `validation/patentes_iba.py` vuelve a comprobar las reglas.

**Consecuencias para el imán (análisis previo a v0.21):**
- El valle calibrado en v0.18 (suelo a 95 mm, Gv = 190) incumple la reivindicación 19. Con el suelo a 100–125 mm (Gv 200–250) el flutter a igual NI sube y ⟨B⟩ baja, así que haría falta todavía más NI.
- Eso aleja aún más el modelo del «< 11 kW» de TUD03T: la tensión potencia/NI empeora.
- Las patentes no la resuelven, pero apuntan a revisar el circuito (corrector de gradiente más profundo, agujeros y retorno de flujo) antes que la curva del acero.

## Revisados en v0.20

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `TUD03T` | Nactergal (IBA), **presentación oral** TUD03 de Cyclotrons 2016 (16 diap.; `c31bd51f…`). Complementa el paper TUD03 | (1) Diap. 4: «Deep valley reinvented (**keeping <11 kW magnet need**)» y 4 patentes (EP16169489 inserto, -90 contorno de polo, -94 corrector de gradiente, -97 geometría). (2) Diap. 5: **mapa Opera de la región central**: dees unidos en el centro por un cuello de ≈ 3–4 mm, puntas a r ≈ 10 mm, 2 chimeneas a r ≈ 12 mm perpendiculares al eje de los dees, movimiento motorizado. (3) Diap. 7: **4 × VHS-6** y bombeo medido (1.9e-5 / 3.6e-6 / < 2.4e-6 mbar a 15/30/45 min tras 2 h abierto). (4) Diap. 8: base 2.7e-7 mbar; transmisión 81 % (150 µA), > 77 % (180 µA), 63 % con 2 ODP; C18/9 55 %. (5) Objetivo > 200 µA; yugo arriba a 1.6 m | **Incorporado.** El vacío se calibra con el bombeo: el acero dulce desgasifica 30 veces más que el metal limpio (antes el modelo era 15–40× optimista); vacío base 2.7e-7. El mapa central **confirma el puente entre dees** de v0.19. Quedan anotados sin aplicar: r_min, chimeneas a 12 mm y la discrepancia de potencia del imán (ver abajo) |
| `KIR30` | Jung et al. (KIRAMS), *Design of the RF system for a 30 MeV cyclotron*, EPAC2006 TUPCH137 (`8d973ccf…`) | Resonador λ/4 con stems verticales, modelo de línea de transmisión + CST, acoplamiento **capacitivo**, sintonizador fino; 64 MHz, Q 6529 | **Validación del método.** Con 2 stems de 99 Ω en paralelo por dee, la línea de transmisión da C_dee = 154 pF frente a 167 pF del FEM (8 %) y 884 A por stem (coherente con el lazo) |
| — | Ha et al., *Design of RF pick-up for the cyclotron*, Cyclotrons2016 THP14 (RFT-30; `c5135960…`) | Sonda coaxial en el valle (campo E), S21 ≈ −75 dB, efecto en f de ±0.03 MHz; las cavidades se diseñan **sobreacopladas sin haz** y críticas con haz | **Referencia.** Confirma el criterio del acoplador de v0.19 (β = 1.17 sin haz). Sonda no modelada |
| — | Brandenburg (KVI), *Cyclotron RF systems*, CAS 2015 (24 diap.; `e058f8e5…`) | Curso: línea de transmisión vs 3D, métodos de sintonía («caja hacia fuera: menos inductancia → f sube», como el pistón de v0.19), acoplamiento capacitivo e inductivo (ventajas e inconvenientes), modos 0/π del PET de GE | **Referencia** (sin datos del KIUBE) |
| — | Medeiros Romão et al. (IBA), *IBA C70 cyclotron development*, Cyclotrons 2007 (`804747d6…`) | RF de IBA: 2 dees y 4 gaps, 30.4 MHz, 65 kV con ≈ 20 kW, amplificador de 100 kW (tetrodo) y modelo CST; 4 criobombas de 5800 l/s | **Referencia de orden de magnitud** (máquina de 70 MeV) |
| — | Jongen et al. (IBA/JINR), *RF system of the cyclotron C400*, Cyclotrons 2007 (`6ba15d73…`) | Cavidades delta de doble gap a 75 MHz, CST sin pérdidas, Q 12 500, ajuste de f escalando los stems; amplificadores de tetrodo de 100 kW | **Desechado** (otra física); solo confirma la práctica de IBA (tetrodos, CST) |

**Discrepancia abierta: potencia del imán.** TUD03T (IBA 2016) dice «< 11 kW magnet need»; IPAC26 da 19 kW. Con 2 t de Cu y la bobina modelada, 11 kW permiten unos 87–91 kA·v, y el FE calibrado (v0.18) pide 116 kA·v = 19.4 kW de bobinas. Este FE es coherente con IPAC26 si los 19 kW son todo bobina, e incompatible con TUD03T salvo que el circuito magnético real pierda mucho menos en el hierro (J_s más alto, tapa o polos distintos). Las patentes EP16169489/-90/-94/-97 (con dibujos) podrían aclararlo.

## Revisados en v0.19

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `VE24` | do Carmo et al. (ICNAS + IBA), *Improved production of novel radioisotopes with custom energy Cyclone Kiube*, Instruments 8 (2024) 38 (10 págs.; `7c874662…`) | **Es el mismo documento que VE24**, incorporado desde v0.6 (planta, trayectorias a 13/18 MeV, horquilla larga, brida en cuña, ventanas de Nb de 75/250 µm, foto de la máquina, colores). Datos que faltaban: cuñas definidas para 13, 14, 15 y 18 MeV; 90 µA en un blanco líquido a 13 MeV; 180 µA en total en doble haz | **Ya incorporado.** Se añade `extraction.ve_energias_cuna`, y `kiube/beam.py` avisa si se pide otra energía |
| `CY11` | Nuttens et al. (IBA), *Design of IBA Cyclone 11 cyclotron magnet*, Cyclotrons 2010 MOPCP072 (3 págs.; `a5713213…`) | Imán de IBA con Opera-3D/TOSCA (simetría de 90°, como nuestro FE). +160 G por 10 A a 186 A; relleno de valle de 50/100 mm → +180/+450 G, con flutter < 5 %. **Curva B-H del proveedor del hierro.** Carrusel de diseño C18/9: horquilla de 25 mm, 2 láminas, 8 azimuts. 42.877 MHz con h = 2 | **Incorporado como validación:** la sensibilidad d ln⟨B⟩/d ln NI vale 0.21 en el C11 y 0.24 en nuestro FE calibrado, lo que confirma un nivel de saturación del hierro parecido. Nuevos `extraction.fork_len` (25) y `foils_por_carrusel` (2) |

## Calibración del imán (v0.18)

No se incorporan documentos nuevos. El circuito magnético se recalibra con las fuentes que ya estaban (`KEY22` Fig. 3, `S09`, `IPAC26`, `T04F1`; detalle en `fem/README_FEM.md`):
- **Acero:** `bh_acero_bajoC_kiube`, con J_s = 2.12 T (CALC).
- **Tapa A:** 181 mm (antes 166, T04F1 ±15).
- **Valle:** 95 mm (antes 110).
- **Tapón:** −16 mm.
- **NI_FE_calibrado:** 116 kA·v, frente a los 170 k anteriores.

## Revisados en v0.17

| Etiqueta | Documento | Qué aporta | Decisión |
|---|---|---|---|
| `AUT15` | *wep09*, ICAP 2015: diseño del imán de un ciclotrón de 18 MeV (AUT) con TOSCA | 4 sectores, colina 46°/valle 44°, polo R 63 cm, gap 3.66–5.99 cm, NI 28 kA·v, acero AISI 1010 (B_sat 1.85 T), B_max 1.72 T. Es una máquina de valle poco profundo | **Solo como referencia:** confirma que el AISI 1010 es la elección típica. Nada más aplicable |
| `BATAN13` | Atom Indonesia 2014: imán del ciclotrón de 13 MeV de BATAN con Opera-3D/TOSCA | NI 88 024 A·v (2 × 44 012), gap de valle 120 mm, B_valle 0.769 T, gaps de colina 40–50.8 mm, 17.3 t, curva B-H de acero de bajo carbono | **Usado como cota, no como objetivo** (corregido en v0.18). Su B_valle/(μ0·NI/g_valle) es 0.84; en el KIUBE calibrado sale 0.58, porque el valle profundo recibe el campo por el borde de las colinas. En v0.17 se dedujeron de aquí unos 100 kA·v: era una extrapolación indebida. La calibración real (v0.18) se hizo con KEY22 y S09 y da 116 kA·v |
| `KEY22` | *Development of the Cyclone KEY* (ya registrado) | Modelo CST; resuena a 41 MHz con 32 kV; simetría perfecta en el plano medio; sin corriente RF en los polos ni liners; valles taladrados para vacío y RF; la cavidad RF es la cámara de vacío y une los 4 agujeros de valle; turbo HiPace 2300 en la cavidad; la Fig. 5 muestra stems altos en una cavidad cilíndrica superior con barra horizontal | **Confirma la topología RF del modelo.** Sin cambio numérico. Indicio de stems acoplados por la barra, anotado en el pendiente RF |

Cada cota de `params.yaml` lleva la etiqueta de su fuente. Este archivo registra los documentos del usuario: qué se tomó de cada uno, cómo de confidencial es y en qué discrepa de lo que ya teníamos.

## Revisados en v0.16

| Etiqueta | Documento | Público | Qué aporta | Decisión |
|---|---|---|---|---|
| `JUAS22` | E. van der Kraaij (IBA), *Accelerators for Medicine and Industry*, JUAS 2022 (174 diapositivas; SHA-256 `f84e4c79…97aef`) | Sí | p.68/157: tabla del Cyclone 18P (= KIUBE): altura 1.85 m, 1.95 × 1.95 m, plano medio 1.12 m, búnker mínimo 3.7 × 3.7 × 2.5 m, hueco de izado 2.2 × 2.2 m, 17 t, hasta 30 Ci, 300 µA. p.44–45: fotos del C18/9 de valle profundo (dees de cobre con **dummy dees a tierra** en los bordes de las colinas, 4 gaps aceleradores, puller a V_dee, **gap chimenea–puller ≈ 1.5 mm**, tapón central, bordes de polo desmontables) | **Incorporado** |
| — | E. Pearson et al. (IBA), *The new IBA superconducting synchrocyclotron (S2C2): from modeling to reality* (5 págs.; `f1aaae41…5f0d529`) | Sí | Sincrociclotrón superconductor de protonterapia (otra física, otra máquina) | **Desechado**: no resuelve ningún pendiente del KIUBE |
| — | S. Gammino et al., *A new superconducting ECR ion source for the future accelerator facilities*, EPAC 2002 (3 págs.; `19c9fda2…bb7e42f`) | Sí | Fuente ECR superconductora de 28 GHz para iones pesados | **Desechado**: el KIUBE usa fuentes PIG internas de H⁻ |

Cambios en el modelo a raíz de JUAS22:
- **Dummy dees** de cobre a tierra en los bordes de las colinas (`dee.dummy_gap` = 6 mm, EST). La RF se recalibró: la caja del dee baja de 97.6 a 86 mm (24 mm bajo el suelo del valle) y la frecuencia queda en 40.680 MHz.
- La foto también es compatible con dees planos. Con placas y dummy dees a 3 mm la frecuencia sale 46 MHz; esta alternativa queda registrada en `out/fem_rf/barrido.json`.
- **Región central del C18/9:** gap chimenea–puller de 2.6 → 1.5 mm. En el KIUBE se mantiene 2.6 mm (EST) y queda anotado el valor del C18/9.
- **Datos ya conocidos, ahora con respaldo público:**
  - plano medio a 1.12 m y altura de 1.85 m, que la versión pública usa en lugar de los del SRG;
  - masa de 17 t (SRG: 18 t; el modelo: 17.5 t);
  - envolvente de 1.95 m.

## Incorporados en v0.12

| Etiqueta | Documento | Fecha | Confidencialidad | SHA-256 (PDF) |
|---|---|---|---|---|
| `SRG` | IBA *Cyclone® KIUBE Site Requirement Guide*, MID 55672 Rev. G, con los planos SRG-0001…0420 (104 págs.) | 01-06-2023 | **Propietario de IBA.** Solo se puede usar para especificar, instalar y mantener el equipo según el contrato; prohibido reproducirlo | `b899b503…6e461ffd` |
| `CBR` | Folleto IBA *Cyclone® KIUBE – Max potential, max capacity* (7 págs.) | jul-2017 | Público (folleto comercial) | `43697cd2…9a66d9` |
| `PROP` | Cotización presupuestaria IBA, Ref. 2025/8499 v02.01, *KIUBE IntegraLab Zumaya Tabasco* (36 págs.) | mar-2026 | **Confidencial** («proprietary and confidential»). Solo se toma la configuración técnica: sin precios, contactos ni condiciones | `d1b26ebf…6a3cc` |

> **Nota (actualizada en v0.17).** Uso personal autorizado por el propietario de los documentos; lo que sigue es solo histórico. Hasta v0.11 el modelo salía solo de literatura pública. Desde v0.12, `params.yaml` contiene cifras del SRG y la configuración del proyecto (PROP). No publiques el modelo ni sus salidas (`params.yaml`, `validation/validacion_srg.json`) sin revisar esto. Los valores `SRG`/`PROP` se pueden sustituir por los anteriores (`DS`, `TUD03`, `EST`) para tener una versión pública.

Los PDF **no** se copian al proyecto, solo se registra su huella. Para volver a verificarlos, colócalos en una carpeta local fuera del repositorio.

### Qué se tomó de cada documento

**SRG (instalación y máquina)**
- **Geometría.**
  - Plano medio a 1120 mm del piso, por lo que `floor_to_top` pasa de 1600 a 1550.
  - Altura total de 1850 mm hasta la cima del cáncamo central M48, que es el único punto de izado superior.
  - Envolvente de 2000 × 2000 mm (máximo de maniobra 2050 × 2050 × 1850).
  - Pies a ±525 mm, con pilares de 220 × 220 y placas de 350 × 350.
  - Masa de 18 000 kg.
- **Autoblindaje** (solo KIUBE 100/150).
  - Mitades de 4720 × 2735 × 2895 mm.
  - Casco de acero de 6.56 t con 34 t de magnetita y 7.5 t de agua, 48.1 t llena.
  - Sala de 6 × 7 × 3.1 m.
- **Búnker.**
  - Interior de 3.8 × 3.8 m como mínimo (recomendado 4 × 4) y 2.5 m de alto.
  - Muros de 2.0 a 2.4 m y techo de 1.7 a 2.3 m, con hormigón de al menos 2350 kg/m³.
  - Foso de 400 mm bajo la máquina.
- **Servicios** (bloque `sitio`): agua enfriada, potencia, calor, aire comprimido, gases, cable RF, ventilación, armarios y distancias máximas.
- **Término fuente y referencias de dosis** (bloque `fuente_rad`):
  - 6.9·10⁻³ n/p y 8.3·10⁻³ γ/p;
  - dosis a 1 m y tras 200/220/230 cm de hormigón (Tablas 8 y 9). Sirven para validar el Monte Carlo.

**PROP (configuración del proyecto)**
- KIUBE 180 (180 µA), 2 fuentes de iones y doble haz.
- 2 blancos Nirta HCFT F18-700 (3.6 ml), 1 blanco HCFT 13N-NH₃ (2.9 ml) y 2 puertas blindadas.
- Se instalan 4 strippers y 4 válvulas de vacío de serie. La asignación a los puertos 1, 3, 5 y 7 es una estimación (EST).
- Superficie del IntegraLab ONE: 440 + 120 m².

**CBR (folleto)**
- Volúmenes de los blancos: inserto de Ga-68 de 3.7 ml y 11C-CO₂ de 55 ml.
- Producción por modelo, por ejemplo KIUBE 180 → 20 Ci de 18F en 2 h. La cotización dice 24 Ci.
- Envolvente del autoblindaje: 5.2 × 4.3 × 2.45 m.

### Discrepancias abiertas

| Dato | Fuente A | Fuente B | Se usa |
|---|---|---|---|
| Rendimiento de neutrones 18O(p,n) | SRG 6.9·10⁻³ n/p | El mismo SRG da 3.55·10¹¹ n/s a 40 µA, es decir 1.42·10⁻³ n/p (incoherente consigo mismo); BON21 (FLUKA) da 4.93·10⁻³ | BON21 por defecto. `normalizacion="fuente_rad"` está disponible en `kiube.mc.target_sources` |
| Mitad del autoblindaje | SRG 4720 × 2735 × 2895 | Ficha DS 4645 × 2770 × 2575; folleto CBR 5.2 × 4.3 × 2.45 m (conjunto) | SRG (el más reciente) |
| Dimensiones totales | SRG 2000 × 2000 × 1850 | Ficha DS y folleto CBR 1.9 × 1.9 × 1.8 m | SRG |
| 18F por 2 h con KIUBE 180 | PROP 24 Ci | CBR 20 Ci | Solo informativo |
| Techo del búnker del proyecto | Plano RPOS (sin altura de techo) | SRG 1700–2300 mm | 1700 (mínimo del SRG) hasta tener el plano |

## Incorporados antes

Ver la cabecera de `params.yaml`: TUP03, TUP04 y sus figuras, TUD03, S09, BON21, PAT, DS, DSR, VE24, KEY22, IPAC26, CERN-CH, CERN-IS, IKON, F7 y RPOS.
