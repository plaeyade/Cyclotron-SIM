# Cotas estimadas (EST) — revisión v0.26 (actualizada en v0.27)

Quedan **36 cotas EST** de 216. Están ordenadas por impacto en la física del modelo (imán, RF, haz, vacío, blindaje). La última columna dice qué documento o medida las cerraría.

**Cerradas o reclasificadas en esta revisión (v0.26):**
- `dee.angle` 38 → **39°** (TUD03T diap. 9): render en planta rectificado con los 4 insertos de polo; 39.1 ± 0.6° promediando los 2 dees. Con él la RF se recalibra: caja del dee a 76.7 mm, 40.656 MHz.
- `vacuum.gas_mass_g_mol` y `temperature_K` → `CONST` (constantes físicas); `vacuum.active_source` y `P0_ultimate_mbar` → `MODELO` (opciones del modelo, no cotas de la máquina).
- `dee.r_max` (500): el render da 475–550 mm. Es compatible, pero sin precisión para cambiarlo.
- La figura VE24 Fig. 3 (CAD en planta) no sirvió: su anillo oscuro queda a ~455 mm con la escala del agujero de bombeo y no corresponde a la cámara (548 mm).
- Cerradas en v0.18–v0.25 (para referencia): región central (puntas, cuello, fuentes, gap, puller, rendija), valle, tapón, tapa, corrector de gradiente y acero.

| Impacto | Bloque | Cota | Valor | Qué afecta | Cómo cerrarla |
|---|---|---|---|---|---|
| alta | coil | `height` | 190.0 | NI ↔ potencia (19 kW) y ⟨B⟩; reparto de cobre | foto o plano de la bobina con escala; masa de cobre por bobina |
| alta | dee | `aperture` | 28.0 | RF (capacidad del dee), enfoque vertical, pérdidas | plano del dee o foto lateral con escala |
| alta | rf | `factor_Q_real` | 0.65 | tamaño del lazo del acoplador y potencia RF (21.8 kW) | Q medido en la puesta en marcha o potencia directa del amplificador |
| alta | vacuum | `h2_sccm_per_active_source` | 2.0 | presión con la fuente encendida; transmisión | consigna de H₂ de la fuente (log de operación) |
| alta | vacuum | `pump_speed_l_s` | 2000.0 | presión de operación y pérdidas de H⁻ por gas | velocidad efectiva de las VHS-6 tras válvula/baffle (catálogo Agilent + plano de la conexión) |
| alta (blindaje) | vault | `floor` | 500.0 | dosis bajo el búnker; Monte Carlo | plano completo del búnker (pendiente nº 2) |
| media | dee | `dummy_gap` | 6.0 | RF (capacidad) y ganancia por gap | foto o plano de la región de los gaps |
| media | dee | `plate_t` | 8.0 | masa del dee homogeneizado, RF (poco) | plano del dee |
| media | dee | `r_max` | 500.0 | RF y extracción (el dee debe pasar el radio de stripping) | render TUD03T da 475–550 mm; un plano lo fija |
| — | extraction | `window_radius` | 590.0 | **sin uso desde v0.27** (los blancos van en la cara del yugo, VE24) | — |
| media | rf | `acoplador_tipo` | lazo | diseño del acoplador (lazo o capacitivo) | foto del interior de la carcasa RF |
| media | rf | `sintonizador_d` | 40.0 | rango de sintonía (±90 kHz por pistón) | foto/plano del sintonizador |
| media | target | `body_d` | 60.0 | término fuente de neutrones/γ cerca del blanco; blindaje local | plano o ficha del Nirta F18-700 / 13N |
| media | target | `body_len` | 95.0 | idem | idem |
| media | target | `window_d` | 20.0 | idem; densidad de potencia en la ventana | idem |
| media | vacuum | `pump_speed_H2O_l_s` | 1300.0 | pump-down (ya calibrado con TUD03T junto con el acero) | idem |
| baja | dee | `punta_t` | 4.0 | campo puntas–tapón (4.4 MV/m) | plano de la región central |
| baja | dee | `r_min` | 60.0 | solo el arranque de la caja; las puntas vienen del mapa de IBA | plano del dee |
| baja | extraction | `arm_h` | 8.0 | cosmético | idem |
| baja | extraction | `arm_w` | 10.0 | cosmético | idem |
| baja | extraction | `carousel_radius` | 28.0 | cosmético/mecánico | foto del carrusel (VE24 Fig. 5 sin escala) |
| baja | extraction | `foil_size` | 10.0 | cosmético | idem |
| baja | ion_source | `inner_len_r` | 160.0 | geometría del vástago de la fuente | plano de la fuente |
| baja | ion_source | `stem_radius` | 16.0 | idem | idem |
| baja | pole | `step2_chamfer` | 8.0 | forma fina del borde de colina | plano del polo |
| baja | rf | `acoplador_az` | 0.0 | orientación del lazo | idem |
| baja | rf | `coax_feed_d` | [40.0, 92.0] | pasamuros 50 Ω (cosmético) | catálogo del pasamuros |
| baja | rf | `sintonizador_off` | 75.0 | idem | idem |
| baja | target | `cavity_d` | 20.0 | volumen ya fijado por la ficha (3.6/2.9 ml) | idem |
| baja | target | `collimator_len` | 40.0 | pérdidas en el colimador | idem |
| baja | target | `collimator_od` | 36.0 | idem | idem |
| baja | target | `puerta_wh` | [260, 260] | sección de la puerta de HDPE (el espesor de 18 cm es dato, VUB24) | ficha de la Nirta Shielded Door |
| baja | extraction | `tubo_haz_d` | [40, 50] | tubo de vacío del puerto | plano de la cámara |
| baja | target | `he_gap` | 3.0 | HCFT no usa He (PROP): solo separación entre láminas | idem |
| baja | vacuum | `seal_exposed_cm2` | 130.0 | idem | idem |
| baja | vacuum | `seal_length_m` | 12.8 | desgasificación de juntas | plano de juntas |
| baja | yoke | `is_port_h` | 80.0 | abertura de la fuente | idem |
| baja | yoke | `port_height` | 120.0 | escotaduras de los puertos (masa de hierro, poco campo) | foto con escala de la cara del yugo |

**Los cinco datos que más cerrarían:**
1. **Plano del búnker** (piso, muros, laberinto, techo): blindaje.
2. **Ficha o plano del blanco Nirta F18-700 / 13N**: término fuente y blindaje local.
3. **Q medido o potencia directa del amplificador**: factor de Q real y acoplador.
4. **Foto con escala o plano de la bobina**: altura de la bobina, y con ella NI frente a potencia.
5. **Velocidad efectiva de las VHS-6 y consigna de H₂**: presión de operación.

Archivo generado al revisar `params.yaml`; las fuentes de cada cota están en sus notas.

**v0.27:** `tgt_window_um` 50 → 35 µm (CON13); puertas blindadas de HDPE de 18 cm (VUB24); blancos en la cara del yugo (VE24 Fig. 3). Nuevas EST de bajo impacto: `puerta_wh` y `tubo_haz_d`.
