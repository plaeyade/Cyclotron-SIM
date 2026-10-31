# Benchmark de blindaje frente al Site Requirement Guide de IBA

Compara el Monte Carlo propio (Geant4) con las Tablas 8 y 9 del SRG (MID 55672 Rev. G).

| Paso | Comando |
|---|---|
| 1. Geometría y configuración | `python mc_benchmark/generar_srg.py` → `out/mc_srg/` (GDML del ciclotrón, sin búnker) |
| 2. Compilar | `c++ -O2 -std=c++17 $(geant4-config --cflags) srg_bench.cc -o srg_bench $(geant4-config --libs) -lG4gdml -lxerces-c` |
| 3. Correr | `mc_benchmark/correr_srg.sh 6000 3.0 "100 150 200 220 230"` (unos 75 min en 2 núcleos) |
| 4. Analizar | `python validation/benchmark_srg.py` → `validation/benchmark_srg.{json,png}` |

Geant4 11.4.3 se instala desde conda-forge con los datos G4NDL 4.7.1:
`micromamba create -p g4 -c conda-forge "geant4=11.4.3=*noqt*" cxx-compiler`

## Modelo

- **Geometría:** ciclotrón completo (LOD neutrónico, teselado conforme) dentro del búnker estándar del SRG: interior de 4 × 4 × 2.5 m, hormigón Portland de 2.35 g/cm³, muros y techo de espesor *t*, losa de piso de 1 m.
- **Fuente:** doble haz con los blancos F18 T1 y T5, 2 × 100 µA (igual que las tablas del SRG). Espectro de 18O(p,n) de BON21 (FLUKA), emisión isótropa.
- **Normalización:** 4.93·10⁻³ n/p (BON21) o 6.9·10⁻³ n/p (SRG).
- **Física:** QGSP_BIC_HP. Las gammas de captura e inelásticas salen del transporte.
- **Dosis:** H*(10) con los coeficientes de ICRP 74, por longitud de traza en volúmenes de aire.
- **Reducción de varianza:** importancia geométrica en capas de 10 cm del muro frontal, factor 3 por capa (mundo paralelo). En los otros muros, la partícula se elimina a más de 40 cm de profundidad.

## Resultados (fuente SRG)

v0.27: los blancos pasan a la cara del yugo (VE24, VUB24), así que la fuente queda a r = 969 mm en lugar de 597 mm. El búnker del SRG se corre **sin puertas de HDPE**, como las Tablas 8–9 (`generar_srg.py --puertas` las añade en `out/mc_srg_puertas/`).

| Punto | v0.26 (r 597) Modelo/SRG | v0.27 (r 969) Modelo [µSv/h] | SRG [µSv/h] | v0.27 Modelo/SRG |
|---|---|---|---|---|
| Frente al blanco (0°), junto al muro | 0.23 | 6.2·10⁷ ± 6 % | 1.2·10⁸ (a 1 m) | **0.51** |
| 90° | 0.54 | 3.3·10⁷ ± 8 % | 2.7·10⁷ | **1.22** |
| 180° | 0.26 | 7.4·10⁷ ± 6 % | 1.2·10⁸ | **0.61** |
| Tras 200 cm | 0.02 | 0.28 ± 21 % | 3.62 | 0.08 |
| Tras 220 cm | 0.008 | 3.1 ± 98 % (un evento gamma) | 0.79 | — |
| Tras 230 cm (1 semilla) | 0.007 | 0.024 ± 32 % | 0.40 | 0.06 |
| TVL corregido por 1/r², ajuste ponderado 100–230 cm | 25.6 cm (3 puntos) | 28.6 cm | 33.9 cm | |

El caso de 230 cm (semilla 2) no terminó en el límite de 2 h y se dejó fuera. En el de 220 cm, el 99.5 % de la dosis viene de una sola gamma de peso alto, así que es ruido estadístico; el ajuste del TVL descarta los casos con error ≥ 100 % y pondera por 1/σ.

## Puertas blindadas de HDPE (v0.27)

`python mc_benchmark/generar_srg.py --puertas` y `correr_caso.sh t semilla 6000 3.0 out/mc_srg_puertas` para 100/150/200 cm × 2 semillas; `python validation/benchmark_puertas.py` → `validation/benchmark_puertas.{json,png}`.

| Magnitud | Sin / con puertas (modelo) | Referencia |
|---|---|---|
| Fluencia de neutrones en el muro frontal, 0–10 cm | **÷1.85–2.0** | IBA: activación a 0–10 cm ÷2.3 (2.7 → 1.15 Bq/g, VUB24) |
| Ídem, 10–100 cm | ÷2.1–2.7 | |
| H*(10) a 0°, junto al muro / a 1 m del blanco | ÷3.3–4.9 | SRG §2.1.8.8: flujo n ÷10 |
| H*(10) a 90° / 180° | ÷2.0–2.6 / ÷2.7–3.3 | |
| H*(10) tras el muro | ÷1.7–9.5 (estadística pobre) | |

Perfil en profundidad con puertas (muro de 200 cm), relativo al tramo 0–15 cm: **1 / 0.45 / 0.13 / 0.035** en 0–15/15–30/30–45/45–60 cm. Las varillas medidas en VUB24 dan 1 / 0.41–0.57 / 0.09–0.19 / 0.07 (Eu-152 y Co-60). El modelo reproduce el factor de IBA para la activación y la caída con la profundidad medida en un KIUBE real. El ÷10 del SRG no sale: probablemente se refiere al flujo hacia adelante justo detrás de la puerta, no a la fluencia en el muro.

## Lectura

1. **Dentro del búnker.**
   - El modelo es coherente con la física básica. A 1 m del blanco da 4.6·10⁷ µSv/h, igual que un punto isótropo analítico (4.7·10⁷).
   - Con la fuente en la cara del yugo (v0.27) el modelo queda a 0.5–1.2 del SRG dentro del búnker (antes 0.23–0.54): **gran parte de la diferencia de v0.26 era la posición del blanco**. El resto (×2 al frente) lo más probable es que IBA modele la emisión de 18O(p,n) **hacia adelante**: el espectro BON21 que usamos es integrado en ángulo y lo emitimos isótropo. Su «a 1 m» también puede referirse a un punto más cercano al blanco.
2. **Tras el muro.**
   - El modelo atenúa más que el SRG: TVL de 28.6 cm frente a 33.9 cm (v0.27). Tras 2–2.3 m queda un factor 12–16 por debajo del SRG (antes 50–140).
   - Nuestra TVL coincide con la teoría de remoción para hormigón ordinario (Σ_R ≈ 0.09 cm⁻¹, TVL de unos 25 cm).
   - La de IBA es la de un hormigón más seco o la de un cálculo conservador.
   - **Sensibilidad:** con la mitad de hidrógeno (0.5 % en peso, hormigón seco envejecido) la dosis tras 150 cm sube unas 10 veces. El contenido de agua del hormigón real es la mayor incertidumbre.
3. **Estadística.** Las dosis tras 2 m tienen un error del 40–60 %, así que valen para el orden de magnitud y no para la cifra fina. Para bajarlo hacen falta más eventos o ventanas de peso (weight windows).

**Conclusión para el proyecto:** para el diseño del blindaje mandan las tablas del SRG y el Monte Carlo de IBA (opción 002.010 de la cotización), porque son conservadores. Este modelo sirve para estudios relativos (orientación de blancos, puertas blindadas, laberinto) y para entender la sensibilidad. No sirve para sustituir el cálculo de IBA ante la autoridad.

## Verificaciones

- **Campo libre a 1 m:** el modelo da 4.64·10⁷ µSv/h y el cálculo analítico 4.68·10⁷ µSv/h (un blanco).
- **Análogo frente a importancia** (muro de 60 cm): la traza por capa coincide dentro de la estadística, así que la división y la ruleta no sesgan el resultado.
- **Balance:** el rendimiento integrado del espectro reproduce 4.93·10⁻³ n/p (BON21), con energía media de 2.61 MeV.

## Techo del búnker del proyecto (v0.28)

| Paso | Comando |
|---|---|
| Configuración | `python mc_benchmark/generar_techo.py --altura 2500` → `out/mc_techo/h2500/` (con puertas; `--sin_puertas`) |
| Compilar | `c++ -O2 -std=c++17 $(geant4-config --cflags) bunker_techo.cc -o bunker_techo $(geant4-config --libs) -lG4gdml -lxerces-c` |
| Correr | `mc_benchmark/correr_techo.sh out/mc_techo/h2500 240 1` (techo grueso: curva por capas; ~1 h) y 150/200 cm para comprobar la salida |
| Analizar | `python validation/techo_bunker.py` → `validation/techo_bunker.{json,png}` |

El techo es un volumen aparte con importancias en capas de 10 cm (+z). La H*(10)-fluencia por capa dentro de un techo grueso es una cota superior de la dosis de salida, porque incluye el retroceso del hormigón de encima. Las corridas reales dan 0.2–0.5 veces la curva, así que el diseño usa la curva sin reducir.
