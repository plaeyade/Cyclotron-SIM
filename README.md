# Cyclotron-SIM

Modelos CAD y físicos de ciclotrones médicos, independientes del fabricante: geometría paramétrica, haz y extracción, vacío (Molflow+), FEM (estructural, magnético, RF) y Monte Carlo (Geant4/OpenMC) del término fuente y del búnker.

> **Uso interno, repositorio privado.** Los documentos fuente son propiedad del usuario, que tiene los derechos de uso.

## Máquinas

| Carpeta | Máquina | Estado |
|---|---|---|
| `machines/iba_kiube` | IBA Cyclone® KIUBE 18 MeV (100–300 µA, doble haz, puerto VE 13–17 MeV) | Modelo v0.36 completo; MC de búnker validado contra la SRG Rev G (F90) |
| `machines/ge_pettrace` | GE PETtrace | pendiente de fuentes |
| `machines/best_cyclotron_systems` | Best Cyclotron Systems | pendiente de fuentes |
| `machines/abt_bg75` | ABT BG-75 | pendiente de fuentes |
| `machines/lbt_longevous` | Longevous Beam Technologies | pendiente de fuentes |

Todas siguen `machines/_schema/maquina.schema.yaml`: cada cota con su fuente o `EST` / `pendiente_fuente`.

## KIUBE: uso rápido

```bash
cd machines/iba_kiube
python -m pytest                       # pruebas del modelo
python build.py --lod mc               # GDML Geant4 + DAGMC
mc_benchmark/compilar.sh               # compila los benchmarks Geant4
```

Ver `machines/iba_kiube/README.md` para todas las salidas (neutronics, fem_struct, fem_mag, vacuum, print, render).

## Monte Carlo del búnker (F89 → F90)

- **F89** incluía en los blancos T1/T5 las *puertas blindadas* opcionales de IBA (blindaje local de HDPE sobre el puerto). Subestimaba la dosis exterior ~10× frente a la SRG.
- **F90** las quita (no se comprarán) y usa doble haz 2×150 µA con el rendimiento IBA de 6.9×10⁻³ n/p. Concuerda con la SRG Rev G: tabla 8 (1 m, 0.79×) y tabla 9 (contacto a 200–230 cm de concreto, IBA/MC 0.5–0.7).
- Configuraciones en `machines/iba_kiube/resultados/mc_f90/`.

## Relación con otros repos

- **radiopharma-os**: gemelo digital de producción (EOB/EOS, QC, métricas). Consume de aquí el término fuente y los rendimientos por blanco (fases F81–F85, también en `rpos_fases/`).
- **Zumaya-AI**: proyecto de instalación (búnker, blindajes, estructura). Usa este modelo como dependencia.
- **RadShieldMx**: cálculo analítico de blindajes. Debe validarse contra los MC de este repo.
- **GAMMA AI**: plataforma maestra.

## Documentos IBA pendientes

Quedan documentos de IBA por integrar: reporte de blindaje Monte Carlo estándar o personalizado, plano de puertos con blancos asignados y datos de puertas blindadas. Van en `machines/iba_kiube/fuentes/`, y cada dato nuevo se cita en `params.yaml`.
