# Molflow+ en el modelo KIUBE

## Qué aporta

[Molflow+](https://molflow.docs.cern.ch/) es el Monte Carlo de partículas de prueba del CERN para
**flujo molecular en ultra-alto vacío**. En un ciclotrón de H⁻ con fuente interna, como el KIUBE, el vacío
fija la fracción de haz que se pierde por *stripping* en el gas residual. Esas pérdidas son, a su vez,
una fuente de neutrones y de activación de la cámara (Bonvin 2021 supuso ~50 % en el 18/9).

Molflow da el **campo de presión 3D** dentro de la cámara. Con él se reemplaza el supuesto de presión
uniforme por el perfil real, que crece cerca de la fuente PIG en el centro, y las pérdidas pasan a
tener distribución espacial y energética.

## Qué genera `python build.py --lod vacuum`

| Archivo | Contenido |
|---|---|
| `molflow/kiube_vacio_cm.stl` | Superficie cerrada del vacío (cm), triángulos ordenados por grupo. Las normales se reorientan HACIA EL GAS por cáscara cerrada (la cavidad hacia dentro, las 8 láminas de stripper hacia afuera) |
| `molflow/grupo_*.stl` | Cada grupo por separado: bombas, paredes por material, fuentes de H₂, ventanas, láminas |
| `molflow/kiube_molflow_H2O_10h.txt` | MolflowCLI, especie H₂O (18 g/mol): desgasificación de metales a 10 h [CERN] y sticking de bombas para agua |
| `molflow/kiube_molflow_H2_fuente.txt` | MolflowCLI, especie H₂ (2 g/mol): gas de la fuente PIG activa y sticking de bombas para H₂ |
| `molflow/kiube_molflow_grupos.json` | Rango de facetas, área y parámetro de cada grupo |
| `perdidas_h0_impactos.csv` | E, x, y, cos θ de los H⁰ que salen de la espiral y golpean la pared (fuente de pérdidas para Geant4) |
| `tpmc_pumpdown.json`, `pumpdown.png` | TPMC propio (misma física que Molflow) sobre el mismo STL: respuesta por material, P(t) del pump-down, H₂ de la fuente y perfil radial |
| `report.json` | Presión (TPMC, Molflow o analítica), σ calibrada, pérdida por MeV y verificación de *stripping* de Lorentz |

Las superficies se clasifican automáticamente por **el volumen vecino** en la malla conforme. Toda
cara que el vacío comparte con acero, cobre o inoxidable es pared de ese material. Las bombas son
la entrada de las 2 difusoras bajo la máquina y la tapa de las 2 carcasas RF; las 4 bombas suman
1774 cm². Las bocas residuales (~37 cm²) quedan cerradas.

## Flujo con Molflow 2.11 (integrado)

`kiube/molflow.py` escribe directamente el XML de Molflow 2.11.1. Ya no hace falta importar el STL en la GUI.

- `molflow/kiube_vacio.zip` es la geometría lista (72 134 facetas, mismo orden que el STL). Se abre en la GUI o en molflowCLI.
- Los archivos `--setParamsByFile` valen tal cual.

```bash
# todo automático: genera la geometría, corre molflowCLI (H₂ de la fuente) y usa SU perfil de presión
# para las pérdidas de H⁻ (busca molflowCLI en MOLFLOW_CLI, en el PATH o en <proyecto>/molflow/)
python build.py --lod vacuum --molflow auto

# o a mano
molflowCLI -f kiube_vacio.zip --setParamsByFile kiube_molflow_H2O_10h.txt -d 1e6 --writeFacetQuantities --outputPath res_H2O
molflowCLI -f kiube_vacio.zip --setParamsByFile kiube_molflow_H2_fuente.txt -d 1e6 --writeFacetQuantities --outputPath res_H2
python build.py --lod vacuum --molflow res_H2/facet_physics.csv
```

**Unidades del XML de Molflow 2.11** (verificadas, ver `kiube/molflow.py`):

- coordenadas en cm;
- `Outgassing constValue` en Pa·m³/s por faceta (1 Pa·m³/s = 10 mbar·l/s);
- `desType = 2` es coseno (`1` es uniforme y baja la transmisión de Clausing un 16 %).

En `--setParamsByFile`, `outgassing` va en mbar·l/s y `specificOutgassing` en mbar·l/s/cm².

En modo dependiente del tiempo, el atributo `molecules_per_TP` de los momentos no está en moléculas por partícula. El lector usa las moléculas reales inyectadas (∫Q dt / kT).

### Molflow 2.11.1 frente al TPMC propio (`validation/molflow_vs_tpmc.py`, figura `validation/molflow_vs_tpmc.png`)

| Caso | Molflow 2.11.1 | TPMC propio | Referencia |
|---|---|---|---|
| Transmisión de tubo, L/R = 1 / 2 / 10 | 0.6722 / 0.5140 / 0.1912 | 0.6704 / 0.5153 / 0.1911 | Clausing 0.6720 / 0.5136 / 0.1910 (1σ ≈ 0.0009) |
| Pump-down del tutorial (meseta, τ) | 4.943·10⁻⁵ mbar, 1.568 ms | 4.936·10⁻⁵ mbar, 1.558 ms | 5·10⁻⁵ mbar; V/S = 1.571 ms |
| KIUBE, H₂O de metales a 10 h (zona del haz) | 9.53·10⁻⁹ mbar | 9.48·10⁻⁹ mbar | +0.5 % |
| KIUBE, H₂ de la fuente, 2 sccm (zona del haz) | 5.69·10⁻⁶ mbar | 5.74·10⁻⁶ mbar | −1.0 % |

Los dos motores coinciden dentro del ruido estadístico.

- **Molflow:** referencia oficial del CERN, la GUI para inspeccionar mapas de presión y 2.1 M impactos/s en 2 hilos. En el KIUBE tarda ~4–5 min por especie con 3·10⁵ partículas.
- **TPMC:** corre sin instalar nada y da la respuesta por material para armar el pump-down P(t) en una sola corrida.

## TPMC integrado (sin Molflow instalado)

`kiube/tpmc.py` implementa el mismo Monte Carlo que Molflow+: moléculas independientes, emisión y re-emisión con
ley del coseno, sticking por faceta y trazado con Embree sobre `kiube_vacio_cm.stl`. Validación:

| Caso | TPMC | Referencia |
|---|---|---|
| Transmisión de un tubo, L/R = 1 / 2 / 10 | 0.6730 / 0.5119 / 0.1901 | Clausing–Cole 0.6720 / 0.5136 / 0.1910 (≤ 0.5 %, 1σ) |
| Recipiente de inox de 1 m², 100 l/s, 10 h de bombeo | 2.98·10⁻⁸ mbar | Chiggiato §2.2: «del orden de 10⁻⁸ mbar» (Q/S = 2.98·10⁻⁸) |
| **Pump-down del tutorial de Molflow+** (tubo L/R = 1, inyección de 1·10⁻⁴ mbar·l/s entre 10 y 20 ms, 1 l/s por tapa), en modo dependiente del tiempo | meseta 4.94·10⁻⁵ mbar; τ = 1.558 ms | Tutorial «Time-dependent examples» (Kersevan & Ady, CERN, 2021): 5·10⁻⁵ mbar; τ = V/S = 1.571 ms (−1.3 % y −0.8 %; R = 1 cm supuesto) |

El modo dependiente del tiempo (`TPMC.run_tof`) es el de Molflow: cada molécula lleva su tiempo de vuelo y su rapidez se sortea
en cada re-emisión de la Maxwelliana de flujo. Script: `validation/molflow_tutorial_pumpdown.py` (figura y JSON en `validation/`).

El flujo molecular es lineal: se corre una vez por material (q = 1) y el pump-down se arma como
P(t) = Σ q_i(t)·R_i. Eso es exacto aquí, porque V/S ≈ 0.1 s es mucho menor que las horas del pump-down.

## Desgasificación por material (CERN)

| Material | Ley | Fuente |
|---|---|---|
| Inox, cobre OFE, aluminio, láminas | q_H₂O = 3·10⁻⁹ / t[h] mbar·l/s/cm² (α = 1) | Chiggiato, arXiv:2006.07124, Ec. 13 y Tabla 3 (3·10⁻¹⁰ a 10 h) |
| Acero dulce del yugo (polos, sectores) | Misma ley; factor de sensibilidad ×3 y ×10 | CERN-CH no lo cubre: **supuesto** |
| Elastómero expuesto de las juntas de Viton | 1·10⁻⁷ mbar·l/s/cm² a 10 h (H₂O), ∝ t^-0.5 | Chiggiato, *Vacuum technology for ion sources*, CAS 2012, arXiv:1404.0960, Tabla 10; ley t^-0.5 de polímeros: CERN-CH Sec. 5 |
| Permeación de agua por las juntas | 1·10⁻⁷ mbar·l/s por toro de Ø6 cm → 5.3·10⁻⁹ por cm; estacionaria a ~2 meses | CERN-CH §6.2 |

## Validación contra el pump-down

- **Caso CERN:** se reproduce exactamente (tabla de arriba).
- **KIUBE, solo con las cargas de la tabla:** 1.3·10⁻⁸ mbar a 10 h y 6.3·10⁻⁹ a 24 h en la zona del haz; S_eff para H₂O ≈ 3350 l/s (geometría v0.7, valle de 110 mm).
  La base medida de 3·10⁻⁷ mbar (TUD03, **con RF**) se alcanzaría en menos de 1 h de bombeo. Es decir, la base medida
  no la fija la desgasificación térmica: hay ~9·10⁻⁴ mbar·l/s más.
  - Candidatos: desorción por RF y multipactor en dees y liners, H₂ residual de la fuente, retroceso de aceite de las
    difusoras y la posición del manómetro.
  - Ni con el acero ×10 se llega a 3·10⁻⁷ tras un día.
  - Esto reemplaza la antigua calibración ×44, que escondía la diferencia.
- **Fuente de H₂ encendida (2 sccm):** 5.7·10⁻⁶ mbar en la zona del haz (TPMC). El analítico daba 5.3·10⁻⁶, y el Cyclone KEY
  midió 1.2·10⁻⁵ con la fuente encendida, en una máquina distinta con 2 turbos.

## Supuestos declarados (reemplazables en `params.yaml › vacuum`)

| Parámetro | Valor | Fuente |
|---|---|---|
| Presión base | 3·10⁻⁷ mbar con RF | TUD03 (medido) |
| Bombas | 4 de 2000 l/s (H₂) y 1300 l/s (H₂O) | TUD03 da 4 ODP; las velocidades son estimadas |
| Desgasificación | Ver tabla CERN | CERN-CH; el acero dulce es supuesto |
| Carga de H₂ por fuente activa | 2 sccm | Estimada (AMIT usa 3.5 sccm) |
| σ(H⁻ + H₂) | ∝ 1/β², magnitud calibrada a 15 % de pérdida | Transmisión 81 % (TUD03); AMIT 14.6 %. Calibración: σ(1 MeV) ≈ 1.0·10⁻¹⁶ cm² |

**Pendiente para un cálculo absoluto:**
- Tabla σ(E) medida (Nakai 1987 / Tabata-Shirai 2000).
- Velocidad real de las bombas.
- Flujo real de H₂ de la fuente.
