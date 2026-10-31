# Modelo de RF del KIUBE (`--lod fem_rf`)

Calcula los modos propios del resonador: el vacío entre conductores con dees, stems, carcasas RF, polos y tapa. Usa NGSolve con elementos de Nédélec de orden 2.

```bash
python build.py --lod fem_rf                 # modo acelerador nominal → out/fem_rf/{report.json, modo_acelerador.png}
python fem/rf_kiube.py --param dee.top_z=90 dee.stem_radius=35 --tag prueba   # estudios paramétricos
python fem/rf_eigen.py --verificar           # casos analíticos → validation/rf_verificacion.json
```

## Planteamiento

- **Dominio:** la mitad superior del vacío (`fem/rf_dominio.py`): cámara, agujeros de los valles a través de C y A, e interior de las carcasas RF, menos todos los conductores. El plano medio actúa como pared magnética, que es la simetría del modo acelerador.
- **Resonador:** cada dee con su stem forma un coaxial λ/4. El stem está en cortocircuito con la tapa de la carcasa y la caja del dee lo carga con su capacidad.
- **Solver:** Arnoldi con desplazamiento e inversión (`ArnoldiSolver`, `sparsecholesky`); el núcleo de gradientes (λ ≈ 0) se descarta.
- **Pérdidas:** Rs = √(ωμ₀/2σ), con todas las paredes de cobre (σ = 5.8·10⁷ S/m). Es una cota optimista: el acero y las juntas reales pierden más.
- **Tensión del dee:** ∫E·dl en el plano medio, desde el eje del dee hasta el centro de la colina vecina.

## Verificación (solución exacta)

| Caso | FEM | Exacto | Error |
|---|---|---|---|
| Cavidad cilíndrica TM010 (R 0.3 m) | 382.4748 MHz | 382.4752 MHz | −8·10⁻⁷ |
| Q de esa cavidad (cobre) | 22 195.0 | 22 195.1 | −2·10⁻⁶ |
| Coaxial λ/4 (L 0.75 m) | 99.9312 MHz | 99.9308 MHz | 3·10⁻⁶ |

## Resultado y calibración

Con la geometría estimada anterior la frecuencia salía muy lejos de la real:

| Geometría | f del modo acelerador |
|---|---|
| Dee de 2 placas de 8 mm, stem Ø90 (antes, EST) | **79.2 MHz** |
| Caja del dee hasta 10 mm del suelo del valle, stem Ø90 | 47.8 MHz |
| Caja del dee con cara superior a 97.6 mm, stem Ø50 (v0.14) | 40.664 MHz |
| Placas del dee + dummy dees a tierra a 3 mm, stem Ø50 (la forma que sugiere la foto JUAS22) | 46.3 MHz |
| Caja a 86 mm + dummy dees a 6 mm + stem Ø50 (v0.16, valle a 110) | 40.680 MHz |
| Caja a 75 mm + dummy dees a 6 mm + stem Ø50 (v0.18, valle a 95) | 40.70 MHz |
| Caja a 72.6 mm + puente central (36 mm) + sintonizador nominal (v0.19) | 40.683 MHz |
| Caja a 75.6 mm, valle a 100 (v0.21) | 40.626 MHz |
| **Caja a 77.9 mm + puntas y cuello del mapa Opera de IBA (v0.22)** | **40.668 MHz** (objetivo 40.646, T04F8) |

**v0.22:** la placa de 36 mm se sustituye por la geometría del mapa Opera de IBA (TUD03T, diap. 5): puntas de dee de 50° hasta r = 8 mm, de 4 mm de espesor, unidas por un cuello de 3.5 mm. El modo en fase sube 0.6 MHz y la contrafase baja de 125 a 110 MHz (el cuello es más inductivo), pero sigue aislada. La caja del dee pasa a 77.9 mm: 40.668 MHz.

**v0.21:** con el valle a 100 mm (patentes de IBA), la caja del dee sube a 75.6 mm, 24.4 mm bajo el suelo: 40.626 MHz (−0.05 %). Aquí la sensibilidad es de −0.21 MHz/mm.

**v0.18:** la calibración del imán subió el suelo del valle de 110 a 95 mm. La caja del dee se recalibró a 75 mm, 20 mm bajo el suelo (antes 24 mm), con una sensibilidad de −0.34 MHz/mm.

**v0.16:** se añadieron los dummy dees a tierra en los bordes de las colinas, que se ven en las fotos del C18/9 de la clase JUAS 2022 de IBA (p.44–45) y forman los 4 gaps aceleradores. Suben la capacidad, así que la caja del dee baja a 86 mm (24 mm del suelo del valle; 1.7 MV/m a 40 kV).

En v0.14, la frecuencia publicada solo se alcanzaba con dos cambios a la vez:
- un dee en forma de caja que casi llena el valle, a 12.4 mm del suelo del valle (3.2 MV/m a 40 kV, por debajo del límite de Kilpatrick, 8.3 MV/m a 40.7 MHz; corregido en v0.19, antes ponía 6.5);
- stems más delgados (Ø50) que dan más inductancia.

Las dos cotas pasan de EST a CALC en `params.yaml`. **El par no es único:** `out/fem_rf/barrido.json` tiene el mapa f(top_z, radio del stem), y otra combinación (más inductancia en la carcasa, un ángulo de dee distinto) también llegaría a 40.646 MHz. La sensibilidad es de unos −0.7 MHz por mm de altura del dee cerca del valor calibrado; las máquinas reales corrigen esto con sintonizadores.

Valores del modo acelerador calibrado:

| Magnitud | Valor |
|---|---|
| Frecuencia | 40.668 MHz (+0.06 %; v0.22) |
| Q (todo cobre) | 5676 |
| Pérdidas a 40 kV | 12.3 kW en la cavidad, para los 2 dees |
| Convergencia de malla (v0.19) | maxh 0.045/gradación 0.4: 40.713; 0.045/0.55: 40.656; **0.035/0.55: 40.833**; 0.030/0.70: 40.816 MHz (con la caja a 72.0). Incertidumbre ≈ ±0.1 % |

**Coherencia con el SRG:** el armario del amplificador RF disipa unos 12 kW de calor. Con las pérdidas de la cavidad (11.4 kW de cobre ideal, más en la realidad) y la potencia del haz (18 MeV × 180 µA = 3.2 kW), la potencia RF del orden de 15–20 kW es coherente en orden de magnitud. No es una validación fina.

## Acoplamiento entre dees, acoplador y sintonizador (v0.19)

```bash
python build.py --lod fem_rf                       # modo acelerador + acoplador (out/fem_rf/report.json → fem_rf.acoplador)
python fem/rf_barrido_v19.py tag rf.sintonizador_gap=200     # cualquier variante (puente, pistón, caja del dee)
python fem/rf_figura_v19.py                        # validation/rf_acoplamiento.png
```

**Por qué hace falta acoplar los dees.** El SRG (§3.6) describe un solo amplificador (armario PS2) y un solo cable de 50 Ω hasta el ciclotrón. En v0.18 los dos resonadores (dee + stem) eran independientes, con los modos en fase y en contrafase degenerados (0.02 MHz). Así, un sintonizador o un acoplador sobre un dee desintoniza solo ese dee: con el pistón en el dee 1 y sin puente, los modos pasan a "solo dee 0" a 40.7 MHz y "solo dee 1" a 42.9 MHz, y el modo acelerador desaparece.

**Puente central** (`dee.puente_w` = 36 mm, EST). Es una placa de cobre que une las puntas de los dees por encima y por debajo de la región central (z = 14–22 mm, entre las chimeneas y el tapón).

| Configuración | En fase | Contrafase | Q (Cu) | P a 40 kV |
|---|---|---|---|---|
| v0.18, sin puente (caja a 75) | 40.70 | 40.70 | 5906 | 11.4 kW |
| Puente, sin pistón (caja a 75) | 39.71 | 123.1 | 5792 | 12.0 kW |
| **v0.19: puente + pistón nominal, caja a 72.6** | **40.683** | **123.0** | **5658** | **12.0 kW** |

El puente añade capacidad sobre el tapón (−1.0 MHz) y lleva la contrafase a unos 124 MHz: el modo acelerador queda aislado, y un único acoplador y un único sintonizador actúan sobre los dos dees por igual (V1/V0 = 0.99).

**Sintonizador** (`rf.sintonizador_*`, carcasa del dee 1). Es un pistón de cobre de Ø40 que baja desde la tapa por el hueco del stem. El FEM muestra que actúa sobre todo como sintonizador inductivo: deriva corriente del stem en la zona coaxial y sube f. El dominio de RF es la mitad superior con espejo en el plano medio, así que el FEM ve una pareja de pistones; por perturbación, el efecto de un pistón real es la mitad.

| Inserción desde la tapa | Δf por pistón | ΔP por pistón |
|---|---|---|
| 0 (retirado) | 0 | 0 |
| 355 mm (punta a 200 mm del dee) | +150 kHz | +0.41 kW |
| 475 mm (punta a 80 mm) | +190 kHz | +0.55 kW |

La posición nominal es la mitad del rango (210 mm de inserción, punta a la salida de la carcasa): **±90 kHz de sintonía**, unas 4 veces el ancho de banda de la cavidad cargada. Basta para la deriva térmica y las tolerancias. Con Ø60 y la punta a 60 mm el efecto era de +0.5 MHz y +1.3 kW por pistón (+1.07 MHz el par del FEM), excesivo. El barrido se hizo con maxh 0.055 m; las mallas con el pistón a menos de 80 mm del dee no caben en 7 GB. Lo que se compara son diferencias.

**Acoplador de potencia** (`fem/rf_acoplos.py`, carcasa del dee 0). Es un lazo inductivo junto a la tapa, donde está el máximo de corriente del resonador λ/4. El lazo está en un plano que contiene el eje del stem, a 40–110 mm de él. Su flujo se integra con el campo del modo propio normalizado a 40 kV: fem = ω·Φ, y P_ext = fem²/(2·Z0).

| Magnitud | Valor |
|---|---|
| Paredes (Q real = 0.65 × Q de cobre) | 18.5 kW |
| Haz (18 MeV × 180 µA) | 3.24 kW |
| **Potencia directa para 40 kV** | **21.8 kW** |
| Altura del lazo para acoplamiento crítico | **30 mm** (área 70 × 30 mm²) |
| fem en 50 Ω | 1475 V |
| Q0 real / Q_ext / Q cargado | 3678 / 3131 / 1691 |
| Ancho de banda (−3 dB) | 24.1 kHz |
| Sin haz (β = 1.18) | Γ = 0.082: refleja el 0.7 % |

**Verificaciones:**
- La corriente del stem en el cortocircuito que implica el lazo (≈ 950 A) coincide con la de una línea de 99 Ω (Ø50 en Ø260) cargada por la capacidad del dee: V/(Z_c·sen βl) con βl ≈ 25°.
- Con 21.8 kW de salida, los ≈ 12 kW de calor del armario del amplificador (SRG §4.4: 80 % de 15 kW) implican una eficiencia del 64 %, típica de un amplificador de tetrodo.

**Campo máximo (v0.22).** Las puntas de 4 mm (z = 14–18) quedan a 9 mm del tapón rebajado (z = 27): 40 kV / 9 mm ≈ 4.4 MV/m, 0.5 Kp. En v0.19, con el puente de 8 mm a 6 mm del tapón, eran 6.7 MV/m, 0.8 veces el límite de Kilpatrick (E_K = 8.3 MV/m a 40.7 MHz, f = 1.643·E_K²·e^(−8.5/E_K)). Los dees de ciclotrón trabajan hasta 1.5–2 Kp, así que es aceptable; con un puente de 4 mm bajaría a 4 MV/m. La cara de la punta del dee (r = 60) frente al tapón (r = 52) da unos 5 MV/m (0.6 Kp).

**Qué falta:** la geometría real del puente, del sintonizador y del acoplador (EST). Es posible que IBA use un acoplador capacitivo en lugar de inductivo; el KEY tiene penetraciones de acoplador y sintonizador por el yugo. El factor de Q real (0.65) es la incertidumbre principal del tamaño del lazo, que escala como √P: con 0.5, unos 34 mm; con 0.8, unos 27 mm.

## Limitaciones

- **Acoplamiento, acoplador y sintonizador (v0.19):** modelados con geometría EST. El lazo del acoplador no se malla (Ø8 mm); su efecto se calcula con el flujo del modo. El pistón sí está en el dominio.
- **Faltan piezas:** las ventanas y el pasamuros cerámico del acoplador.
- **Simplificaciones:** las fuentes de iones, los strippers y los blancos se omiten. Las cavidades de las bombas inferiores se tratan como simétricas de las bridas ciegas superiores; a 40 MHz los agujeros están muy por debajo del corte, así que su efecto es pequeño.
- **Validación pendiente:** compararlo con un plano real del dee y del stem, o con la frecuencia, Q y potencia medidas en la puesta en marcha.
