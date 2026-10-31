# KIUBE — modelos Monte Carlo (v0.3)

Todo sale de **una misma malla de superficie conforme** generada con Gmsh *fragment*: las caras compartidas entre volúmenes tienen triángulos idénticos. Por eso DAGMC y GDML describen exactamente la misma geometría.

## Archivos

| Archivo | Código | Unidades | Notas |
|---|---|---|---|
| `kiube_openmc.h5m` | OpenMC + DAGMC | cm | 121 volúmenes; `mat:<nombre>`; vacío = `void`; aire = implicit complement |
| `kiube_dagmcnp.h5m` | DAG-MCNP6 (también FluDAG/FLUKA) | cm | + graveyard; `mat:N/rho:-ρ` (IDs en `dagmcnp_materials.json`) |
| `kiube_geant4.gdml` | Geant4 (`G4GDMLParser`), ROOT | mm | 121 `G4TessellatedSolid` + materiales elementales/isotópicos embebidos |
| `sources.json` | todos | mm | 8 fuentes en el centro de cada cavidad + espectro H₂¹⁸O |
| `materials.yaml` | todos | — | composiciones con trazas de activación (Co, Ag, Eu) |
| `run_openmc.py` | OpenMC ≥ 0.14 con DAGMC | | flujo, H*(10), tasas de activación |
| `run_geant4.py` | geant4_pybind ≥ 11.2 | | fluencia por volumen; `--geantino` para prueba sin datasets |
| `validate_g4.py` | geant4_pybind | | solapes del navegador + geantinos |

## Qué se verificó aquí

| Prueba | Resultado |
|---|---|
| Solapes entre los 121 sólidos CAD (intersección booleana) | 0 |
| Volumen teselado / volumen CAD, cada volumen | 0.9936 – 1.0007 |
| Aristas no-manifold u orientación inconsistente | 0 (se corrigió una arista cuchillo ranura–tapón) |
| `G4PVPlacement::CheckOverlaps` (2000 puntos/volumen, tolerancia 1 µm) | 0 de 121 |
| Geantinos desde las 8 cavidades (20 000) | 100 % nacen en el agua enriquecida; 0 atascados; caminos simétricos arriba/abajo |
| `.h5m`: grupos de material, volúmenes, extensión en cm | OK (inspección con h5py) |
| `run_openmc.py` | **no ejecutado** (OpenMC no se instala aquí); sintaxis verificada |

Un hallazgo de esta validación: con 32 segmentos por vuelta, el polígono de la bobina mordía la pared de la cámara (hueco de 2 mm) hasta 0.68 mm. Ahora la malla se refina en los anillos r = 535–585 mm y 745–790 mm, con flecha ≤ 0.3 mm.

## Término fuente

- **Espectro:** H₂¹⁸O (98 %) con protones de 18 MeV, de Bonvin et al. 2021, Fig. 7 (FLUKA). Lo extraje de los vectores del PDF: 238 puntos, de 5·10⁻⁹ a 14.7 MeV.
- **Rendimiento:** 4.93·10⁻³ n/protón. Con 150 µA repartidos en 2 blancos simultáneos, cada uno emite 2.31·10¹² n/s.
- **Corrección empírica:** Bonvin necesitó un factor **×2.5** para reproducir la activación medida en las bobinas del 18/9, atribuido a la emisión hacia atrás subestimada. Se activa con `--bonvin-correction` en `run_openmc.py`.
- **Angular:** la fuente es isótropa. Una distribución angular mejor (Hagiwara 2011, Bakhtiari 2020) reduciría la incertidumbre principal.
- **Pendiente:** las pérdidas de haz en la cámara (~50 % en el 18/9, Bonvin) no están incluidas. Dominan solo donde la activación es baja.

## Uso

```bash
# OpenMC (conda install -c conda-forge openmc=*=dagmc*)
python run_openmc.py --targets T1 T2 --particles 200000 --batches 20

# Geant4
pip install geant4_pybind          # descarga datasets la primera vez
python run_geant4.py --targets T1 T2 --events 100000 --physics QGSP_BIC_HP

# DAG-MCNP: usar kiube_dagmcnp.h5m con `dagmc` en el input; materiales m1…m13 según
# dagmcnp_materials.json (con expansión isotópica para ENDF/B-VIII.0)
```

## Otros códigos

- **FLUKA:** FluDAG lee el mismo `.h5m`. Para FLUKA nativo (CSG) habría que reconstruir la geometría con cuerpos; se puede generar desde `params.yaml`, simplificando polos y blancos.
- **PHITS / MCNP-UM:** usan malla tetraédrica, que cad-to-dagmc puede generar con `unstructured_volumes` sobre volúmenes elegidos (bobinas, yugo) para tallies de activación por elemento.

## Rendimiento esperado

Según Du (MIT 2024), DAGMC con MOAB corre 2–3× más lento que CSG en OpenMC, y con Embree (double-down) queda al nivel de CSG. Este modelo tiene 4.8·10⁵ triángulos, del orden de sus modelos de prueba.
