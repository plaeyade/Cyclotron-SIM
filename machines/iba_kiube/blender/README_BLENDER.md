# Cyclone KIUBE → Blender → impresión 3D

## Flujo en dos comandos

```bash
# 1) CAD → cuerpos imprimibles a escala (STL + .glb + manifest.json con validación)
python build.py --lod print3d --escala 0.1 --impresora fdm        # 1:10, máquina de 174 mm
#   --escala 0.2 → 1:5 (348 mm: el script de Blender divide en cuadrantes)
#   --impresora resina → espesor mínimo 0.3 mm en lugar de 0.8 mm

# 2) Blender headless: pines de registro, orientación, cama, STL por pieza y por placa, .blend
blender -b -P blender/blender_prep_impresion.py -- \
    --manifest out/print3d/manifest.json --salida out/print3d/blender --cama 220x220x250
```

Probado con Blender 5.0.1 (`bpy`) y compatible con 4.2+. Usa el booleano **Manifold** (4.5+) y cae a **Exact** en versiones anteriores.

## Qué produce

| Archivo | Uso |
|---|---|
| `out/print3d/kiube_print.glb` | Abrir en Blender (File › Import › glTF): máquina armada con nombres y colores por material |
| `out/print3d/stl/*.stl` | 41 cuerpos (1:10) en posición de ensamble |
| `out/print3d/manifest.json` | Escala, piezas por cuerpo, estanqueidad, error de volumen, espesor mínimo, orientación |
| `out/print3d/blender/kiube_impresion.blend` | Colecciones «Ensamble» (armado) e «Impresion» (orientado en cama) |
| `out/print3d/blender/piezas/*.stl` | Cada pieza orientada y apoyada en z = 0, lista para el slicer |
| `out/print3d/blender/placas/placa_NN.stl` | Piezas acomodadas por placa de impresión |
| `out/print3d/blender/informe_impresion.json` | Chequeos bmesh, pines, avisos y placas |

## Cuerpos imprimibles

Las 79 piezas del LOD `print` se agrupan en 41 cuerpos:

- **Yugo inferior y superior:** cada uno es una sola pieza, unión de la tapa A, el yugo lateral B, los 4 sectores C, los 4 polos D y el tapón F.
- **Insertos E:** 8, aparte, para imprimirlos en otro color y encajarlos en la ranura.
- **Resto, por separado:** bobinas (2), cámara (1), dees con su stem (4), 8 blancos (cuerpo + colimador), 8 strippers (brazo + carrusel), 2 fuentes, 4 patas y 2 gatos.

## Verificación hecha (1:10, FDM)

- **Estanqueidad:** 41/41 cuerpos estancos y con orientación consistente. Error de volumen frente al CAD ≤ 0.25 %.
- **Después de Blender:** 41/41 STL estancos, comprobado de forma independiente con trimesh, y 0 aristas no-manifold y 0 caras degeneradas (bmesh).
- **Pines de registro:** 4 de Ø4 × 3 mm sobre las diagonales del yugo lateral B. Volumen verificado: +150.6 mm³ macho y −197.2 mm³ hembra, con 0.15 mm de holgura radial.
- **Placas:** 6 en una cama de 220×220 mm. Las mitades del yugo van solas, tapa abajo, con polos y sectores creciendo hacia arriba sin soportes.
- **Escala 1:5:** yugo, bobinas y cámara se dividen en 4 cuadrantes; resultan 56 piezas en 23 placas, todas estancas.

## Avisos a 1:10 con FDM

Las placas de los dees y los brazos de los strippers quedan en **0.8 mm**, justo el mínimo con boquilla de 0.4. Opciones:

- imprimir en resina (`--impresora resina`, mínimo 0.3 mm);
- subir la escala, por ejemplo `--escala 0.125` (1:8);
- engrosarlas en Blender (`Solidify`) aceptando que dejan de ser exactas.

Las astillas locales de 0.07–0.3 mm (bordes de puertos y perforaciones de la cámara) no afectan la pieza: el percentil 2 del espesor queda en 1.75–4.6 mm.

## En la interfaz de Blender

1. Abre `kiube_impresion.blend`. Las unidades ya están en mm (`scale_length = 0.001`).
2. Para revisar, opcionalmente, instala la extensión *3D Print Toolbox* (Edit › Preferences › Get Extensions). Desde 4.2 no viene incluida. Sus chequeos de *Thickness* y *Overhang* complementan el informe.
3. Exporta con File › Export › STL (`.stl`) o con el exportador 3MF de las extensiones, si tu slicer lo prefiere.

> Los pines y los cortes por cuadrante existen solo para imprimir. No forman parte del modelo físico y el informe los declara.

## Render fotorrealista 4K/8K

```bash
python build.py --lod render                       # 200 piezas + out/render/scene.json
# 4K (3840×2160) en una pasada, ~30 min CPU con 96 muestras
blender -b -P blender/render_kiube.py -- --escena out/render/scene.json \
    --salida out/render/kiube_4k.png --res 3840x2160 --muestras 96 [--vista abierta] [--blend escena.blend]
# 8K (7680×4320) con poca RAM: 2×2 regiones con solape y fundido; cada región pesa como un 4K
blender -b -P blender/render_kiube.py -- --escena out/render/scene.json \
    --salida out/render/kiube_8k.png --res 7680x4320 --muestras 48 --tiles 2 --solape 128
python blender/unir_tiles.py out/render/kiube_8k_tiles.json out/render/kiube_8k.png
```

En una pasada, el 8K necesita más de 7 GB de RAM. En este contenedor el OOM-killer lo mató, y por eso existe `--tiles`.

### Con GPU NVIDIA (por ejemplo, RTX 5090 Laptop, 24 GB)

```bash
blender -b -P blender/render_kiube.py -- --escena out/render/scene.json \
    --salida out/render/kiube_8k.png --res 7680x4320 --muestras 128 --gpu auto
```

- `--gpu auto` prueba OptiX, luego CUDA, HIP, oneAPI y Metal. El log dice qué GPU usó; si no encuentra ninguna, sigue en CPU.
- Con 24 GB de VRAM, el 8K entra en una pasada: no hace falta `--tiles`. Puedes subir a 128–256 muestras.
- Las RTX 50 (Blackwell) necesitan Blender 4.4 o posterior (probado con 5.0) y un driver NVIDIA reciente.
- En una laptop, conecta el cargador y usa el perfil de máximo rendimiento. Con batería la GPU baja mucho su potencia.
