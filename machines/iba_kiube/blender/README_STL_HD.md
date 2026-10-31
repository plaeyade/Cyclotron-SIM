# KIUBE en STL de alta definición: por pieza y completo

Se genera con `python build.py --lod stl_hd --escala 0.1 --calidad alta`. Las opciones son `--escala 0.1|0.2|1.0`, `--calidad media|alta|ultra`, `--remallar <mm>` y `--puente <mm>`.

## Qué hay en cada carpeta

| Variante | Carpeta | Triángulos (61 piezas) | Triángulos (completo) | Tamaño [mm] |
|---|---|---|---|---|
| 1:10 alta | `escala_1-10_alta/` | 166,608 | 205,030 | 216 × 192 × 183 |
| 1:10 ultra | `escala_1-10_ultra/` | 525,566 | 618,010 | 216 × 192 × 183 |
| 1:5 alta | `escala_1-5_alta/` | 178,046 | 216,978 | 432 × 384 × 366 |
| 1:1 alta (maestro) | `escala_1-1_alta/` | 286,968 | 323,478 | 2160 × 1920 × 1830 |

- `piezas/*.stl`: los 61 cuerpos imprimibles, que cubren las 196 piezas del modelo. Cada uno es un sólido cerrado, en posición de ensamble (origen en el centro de la máquina, plano medio z = 0).
- `completo/kiube_completo.stl`: **un solo sólido**, la unión booleana real de los 61 cuerpos con manifold3d. Para imprimirlo de una pieza se añadieron *puentes de impresión* (ver abajo).
- `completo/kiube_completo_multicuerpo.stl`: los 61 cuerpos juntos sin unir, sin puentes.
- `completo/kiube_completo.3mf`: un objeto por cuerpo, con los colores de la paleta de las fichas (v0.10). Los cuerpos con varios colores van como objeto con partes, una por color. Es el mejor para PrusaSlicer, Bambu Studio, OrcaSlicer o Cura: asignas un filamento a cada parte. (En v0.9 este 3MF salía sin color: trimesh no exporta colores a 3MF. Ahora lo escribe `kiube/palette.py`, validado con lib3mf en modo estricto.)
- `completo/kiube_completo_color.glb`: un material PBR por color (visores web, Blender, Windows).
- `por_color/`: los cuerpos multicolor partidos en sólidos estancos, uno por color (yugo blanco + polos de acero, fuente Cu + inox, blanco Nb + colimador Al, actuador + tapa).
- `lista_colores.csv`: cuerpo → color (hex, RGB, RAL aproximado, filamento sugerido, fracción de área, confianza).
- `manifest.json`: por cuerpo, triángulos, estanqueidad, variedad (manifold), auto-intersecciones reales, error de volumen frente al CAD, espesor mínimo, si es imprimible con la boquilla elegida, y los puentes.
- Con 1:10 se incluye `blender/` (pines de registro, orientación, 8 placas de 220 × 220 mm). A 1:5 conviene correr `blender/blender_prep_impresion.py`, que parte el yugo, las bobinas y la cámara en cuadrantes para la cama.

## De dónde sale el detalle (y su límite)

- **Origen de los triángulos:** salen directamente de las superficies exactas del CAD paramétrico, con flecha ABSOLUTA a la escala de salida (verificada sobre un cilindro: la flecha real queda en la mitad de la pedida).
  - **alta:** flecha 0.005 mm y ≤ 0.05 rad entre facetas (≥ 126 segmentos por vuelta).
  - **ultra:** flecha 0.002 mm y ≤ 0.02 rad (≥ 314 segmentos por vuelta).
  - Una impresora FDM resuelve ~0.05–0.1 mm y una de resina ~0.02–0.05 mm. Con "alta" las facetas ya no se ven al imprimir.
- **Límite:** más polígonos no agregan detalle que el modelo no tenga. La geometría es la reconstrucción del KIUBE desde literatura pública; las cotas marcadas EST en `params.yaml` siguen siendo estimadas.
- **Las caras planas usan pocos triángulos largos.** Es lo normal en CAD, no afecta al slicer, y por eso el total es moderado.

## Remallado: está implementado, pero para imprimir no hace falta

`--remallar 0.5` aplica el remallado isotrópico de MeshLab (pymeshlab), preservando las aristas vivas. Medido aquí a 1:10:

| Cuerpo | Triángulos | Desviación de Hausdorff | Δ volumen |
|---|---|---|---|
| bobina_sup (arista 0.5 mm) | 2 892 → 297 006 | < 0.0001 mm | 0.000 % |
| dee_0_sup (0.3 mm) | 1 256 → 215 526 | < 0.0001 mm | −0.008 % |
| yugo_superior (1.0 mm) | 10 614 → 258 012 | < 0.0001 mm | +0.006 % |

Multiplica los triángulos por 25–170 sin cambiar la forma: subdivide las facetas que ya existen y no puede recuperar una curva más fina. La suavidad real se obtiene bajando la flecha del CAD (`--calidad ultra`).

El remallado sirve para otra cosa: FEM con elementos regulares, esculpir o pintar texturas. En programas externos equivalen:

- MeshLab: *Remeshing › Isotropic Explicit Remeshing*.
- Blender: modificador *Remesh*, modo *Voxel*. Ojo: el modo vóxel redondea las aristas vivas.
- Instant Meshes: mallas de cuadriláteros.

## Validación (todas las variantes)

- 61/61 cuerpos estancos, con orientación consistente y variedad (manifold3d).
- 0 auto-intersecciones reales. MeshLab marca candidatas en caras planas con triángulos muy finos; una prueba exacta (área de la unión igual a la suma de áreas en cada plano) confirma que no se solapan.
- Error de volumen frente al CAD ≤ 0.011 % (alta) y ≤ 0.002 % (ultra).
- `kiube_completo.stl`: 1 componente, estanco y variedad.
- Corrección en el CAD hecha para esto: en los blancos, el taladro de la ventana y la cavidad (misma Ø) ahora son un solo taladro, sin la cara plana de ancho cero que dejaba OpenCASCADE.

## Imprimir el completo de una pieza

- **Puentes de impresión.** En la máquina real, algunas partes cuelgan de soportes que el modelo no tiene. Para que el completo sea un solo sólido se unen al bloque con puntales de Ø 2.4 mm a 1:10 (`--puente 1.2`; `--puente 0` los quita). Son 15:
  - 2 gatos de elevación con sus 4 auxiliares, separados 0.6–5.9 mm;
  - el anillo de tuberías, a 8.1 mm;
  - los 8 carruseles de stripper, a 0.1 mm.

  No son geometría de la máquina y van declarados en `manifest.json`.
- **Partes interiores.** Dees, polos, cámara y blancos quedan encerrados dentro del yugo: impresos de una pieza no se ven y el slicer los rellena o les pone soporte. Para mostrar el interior, imprime por piezas (`piezas/` o el 3MF), o las dos mitades del yugo por separado.
- **Piezas al límite de la boquilla.** A 1:10 con FDM, las placas de los dees, los brazos de stripper y las difusoras quedan en 0.8 mm (`imprimible: "limite"` en el manifest). Usa resina (`--impresora resina`) o escala 1:5.

## Colores (v0.10)

La paleta sale de las fichas técnicas IBA: `palette.yaml`, generada por `palette_src/muestrear.py`. Las regiones medidas se ven en `docs/paleta_regiones_*.png` y la carta en `docs/paleta_kiube.png`.

- Variante `foto` (por defecto) = máquina real. Variante `render` = render oficial (`--paleta render`). Variante `didactico` = colores por material de v0.9.
- La banda gris de los chaflanes del yugo es color **de cara**, no un sólido aparte. Se ve en el 3MF y el GLB, pero para imprimirla en otro color hay que pintarla en el slicer (pintura multicolor) o a mano.
- Las fichas no dan códigos RAL. El RAL es el más cercano a lo medido y sirve solo de orientación. Compáralo con una carta física antes de pintar una maqueta.
