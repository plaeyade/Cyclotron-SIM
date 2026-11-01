# Gemelo físico KIUBE (fases trasladadas desde RadioPharma OS)

Física de la máquina, del vacío y del búnker que se desarrolló en RadioPharma OS (F63, F65–F79 y
re-corridas v0.35 de F81–F85), trasladada aquí el 2026-11-01. Detalle en `docs/FASES.md` de
RadioPharma OS.

- `pipelines/`: una fase por archivo, con sus pruebas y `conftest.py`.
- `app/domain/`: `gas_stripping`, `vacuum_tpmc`, `orbita_central`. Se conserva el paquete `app`
  para no cambiar los imports.
- `data/kiube/`: perfil ⟨B⟩(r) y flutter medidos (KEY22, Cyclotrons2022).
- `tools/`: generador de cámara Molflow y sellado DAGMC.

Evidencia externa bajo `RPOS_AGENTS_ROOT`. Sin ella, las pruebas que la necesitan se omiten con su
motivo; con `RPOS_REQUIRE_EVIDENCIA=1` se exigen.

```bash
cd machines/iba_kiube && python -m pytest gemelo_fisico/pipelines   # 22 pasan, 36 omitidas sin evidencia
```
