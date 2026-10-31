# Cyclone® KIUBE 18 MeV — modelo CAD paramétrico (v0.36)

Modelo reconstruido a partir de literatura pública y de los documentos aportados. No es la geometría de fabricación de IBA.

> **Uso personal, no público.** Los documentos fuente (SRG, PROP, papers) son propiedad del usuario, que tiene los derechos de uso; el modelo usa toda su información técnica sin restricciones. `herramientas/version_publica.py` queda solo como utilidad opcional.

Cada cota de `params.yaml` lleva su fuente. Las cotas **medidas en figuras** tienen etiqueta propia (`T04F2` = planta a escala de TUP04, ±3 mm; `S09` = mapa de campo medido). Lo **estimado** va como `EST`.

## Uso

```bash
python -m pytest                 # pruebas rápidas (~5 min): fuentes, geometría, masas, vacío, paleta, 3MF, SRG
python -m pytest -m ""           # todas (~8 min), incluye RF contra soluciones exactas y frecuencia del modo acelerador
python herramientas/version_publica.py   # copia sin datos SRG/PROP en dist/kiube_publico (y la prueba)
pip install cadquery pyyaml pyvista
python build.py --lod neutronics   # 135 sólidos (láminas µm, He, H2-18O) → DAGMC
python build.py --lod fem_struct   # 95 sólidos ≥0.5 mm, 1 STEP por pieza → Abaqus/ANSYS/COMSOL SM
python build.py --lod fem_mag      # hierro + bobina + aire conforme (sector 90° y 360°) → COMSOL AC/DC/Opera
python build.py --lod print        # 196 STL (incluye exterior cosmético)
python build.py --lod stl_hd --escala 0.1 --calidad alta   # STL HD: piezas + completo (1 sólido) + 3MF/GLB con color
python palette_src/muestrear.py    # re-mide la paleta de colores en las fichas → palette.yaml (--paleta foto|render|didactico)
python build.py --lod print3d --escala 0.1   # 61 cuerpos imprimibles → Blender (blender/README_BLENDER.md)
blender -b -P blender/blender_prep_impresion.py -- --manifest out/print3d/manifest.json --salida out/print3d/blender
python build.py --lod vacuum       # superficies del vacío → Molflow+ y pérdidas de H⁻ (out/vacuum/README_MOLFLOW.md)
python build.py --lod render       # 196 piezas con material visual → blender/render_kiube.py (4K/8K)
python build.py --lod mc           # 137 volúmenes: .h5m OpenMC / DAG-MCNP + .gdml Geant4 (ver out/mc/README_MC.md)
python render.py && python render3d.py
python validation/validar_documentos.py     # modelo vs VE24 / KEY22 / IPAC26 → validation/validacion_v06.{json,png}
python validation/validar_srg.py            # modelo vs Site Requirement Guide (cotas de instalación, masa) → validation/validacion_srg.json
python validation/figura_ve.py 13           # planta del puerto de energía variable a 13 MeV
```

Energía variable: en `params.yaml` pon `extraction.ve_energy_MeV` entre 13 y 17 y vuelve a construir cualquier LOD. Los puertos 3 y 7 pasan a la horquilla larga, con el blanco girado y la brida en cuña.

## Un modelo, cinco salidas

| LOD | Para | Contenido | Formato |
|---|---|---|---|
| `neutronics` | OpenMC/MCNP vía DAGMC | todo, incl. láminas de 5–50 µm, He, agua enriquecida | STEP por material + `cubit_dagmc.py` |
| `fem_struct` | Abaqus, ANSYS Mechanical, COMSOL Structural/Heat | hardware ≥0.5 mm, sin fluidos ni láminas | STEP por pieza + `parts_materials.json` (E, ν, k, cp, α, σ) |
| `fem_mag` | COMSOL AC/DC, Opera-3D, ANSYS Maxwell | solo hierro + bobinas + **dominio de aire conforme** (la suma de volúmenes = caja, verificado) | STEP por dominio + `domains.json` + `SETUP_FEM_MAG.md` |
| `print` | impresión 3D | igual que fem_struct | STL por pieza |
| `mc` | OpenMC, DAG-MCNP, FluDAG, Geant4 | + vacío explícito y búnker de concreto (EST); malla conforme compartida, validada en Geant4 | `.h5m` ×2, `.gdml`, `sources.json`, scripts |
| `fem_rf` | NGSolve (y STEP para COMSOL RF / CST) | vacío entre conductores (mitad, pared magnética en el plano medio) | `dominio_nominal.step` + `report.json` + figura del modo |

Todos salen del mismo `params.yaml`: cambiar una cota regenera los cinco.

## Qué cambió en v0.36 (blindaje del búnker del proyecto con la v0.35: equivalente de RPOS F70)

- **`mc_benchmark/bunker_mapa`** (Geant4, QGSP_BIC_HP): la máquina v0.35 en el bloque BIND de F70 (sala 5.9 × 4.7 × 2.5 m, muros 2.11 m, losa 0.55 m, techo 1.45 m), con los blancos F18 en la cara del yugo y las puertas de HDPE. Fuente BON21 a 2 × 90 µA.
  - Malla de H*(10) de 10 cm (n + γ, recorrido exacto de vóxeles) y fluencia por región (térmica 1/v, epitérmica, rápida) para la activación.
  - Cascarones laterales de importancia cada 10 cm a través de los muros.
  - Uso: `python mc_benchmark/generar_mapa.py`; luego `mc_benchmark/correr_mapa.sh semilla paso segundos [f]` (bloques de 2000 historias, robustos a historias muy largas); luego `python mc_benchmark/analizar_mapa.py out/mc_f70v035/h2500_t1450`; figura con `mc_benchmark/figura_mapa.py`.
- **Resultados** (12 corridas × 2000 historias, error estándar entre corridas):
  - TVL del muro 25–28 cm (benchmark SRG v0.27: 28.6; SRG Tabla 9 ≈ 30);
  - **dosis exterior directa** de 0.002–0.008 µSv/h (media de la cara, 0–30 cm fuera) y 0.027 ± 0.010 µSv/h frente a T1; transmisión del muro 0.6–2.5·10⁻⁹;
  - fracción de la dosis fuera del bloque 2.7·10⁻⁴;
  - activación de saturación: Co-60 10 GBq, Eu-152 4.8 GBq, Cu-64 590 GBq, Ag-110m 0.57 GBq.
- **Contra RPOS:** F70 no tenía cuentas fuera. F78 las extrapoló con μ = 0.0142 cm⁻¹ (TVL 162 cm), lo que sobreestima la exterior unos 10⁷. Las activaciones de F72 son del mismo orden salvo la de Ag-110m.
- **Limitación:** meseta de varianza a 1.4–2.0 m de profundidad (neutrones que cruzan el techo sin sesgo). El ajuste de μ usa 40–130 cm.
- `validation/f70_v035_mapa.png`. Radiopharma-OS: fase F83 (parche).
- **Cierre de la cadena F70–F78** (Radiopharma-OS F85):
  - 5 corridas más con `mc_benchmark/bunker_mapa_v1b` (mismo esquema de importancia y malla de fluencia de neutrones en 3 grupos): 17 corridas, 34 000 historias;
  - `analizar_mapa.py` exporta `mallas_gemelo_f77.npz` (H*(10) media/σ y fluencia en 3 grupos, formato F77);
  - un esquema v2 con importancia también en el techo y la losa se descartó: no completaba 500 historias en 25 min (`mc_benchmark/bunker_mapa_v2_descartado.cc`).
- **Molflow+ 2.11.1** (paquete Debian, `molflowCLI`): validación del vacío (Radiopharma-OS F84). Con el caso analítico se encontraron y corrigieron tres errores en Radiopharma-OS:
  - el motor TPMC daba presiones ×2;
  - el generador de Molflow dejaba una corona abierta en la tapa, por la que escapaba el gas;
  - la fuente de H₂ inyectaba solo la mitad.

  Con las correcciones, Molflow y el TPMC coinciden al 1 % con la cámara v0.35, y la pérdida por gas es del 4.4 % (2 sccm; 7.0 % con 3.5 sccm).

## Qué cambió en v0.35 (simulación de punta a punta: fuente → blanco)

- **`fem/fuente_a_blanco.py rc | acel | ext`**: las tres etapas ven el mismo campo, el mapa FE oficial × nivel, con f_RF = 40.640 MHz.
  1. **Región central** con el menisco autoconsistente (puller con bisel, chimenea Ø6) y el campo E de la FE; 800 partículas emitidas en −80…+40°. Salen en r = 126 mm (≈1.25 MeV) con posición, momento (ajuste cuadrático de los 3 últimos puntos), energía y fase RF.
  2. **Aceleración** (`fem/aceleracion.acelerar`, ahora con `estado=` y `lamina=`) hasta la lámina del puerto 1: azimut 32.2°, desde la órbita cerrada de 18 MeV en ese azimut (463.0 mm; la óptica de v0.2 daba 463.4) más 10 mm.
  3. **Extracción**: H⁻ → p con la misma velocidad y dispersión múltiple en 5 µm de C (1.1 mrad). Los protones se siguen en el campo de borde real hasta la cara del yugo y luego en recta hasta la ventana del blanco. El azimut de la lámina se reajusta para centrar el haz, como hace el operador con el carrusel.
- **Mapa extendido del imán** (`bz_midplane_ext_NI115k.npz`, 0.40–0.95 m, con el material de cada punto). El mapa anterior quedaba fijo a partir de 0.58 m, porque el muestreo fuera de la malla devolvía el último valor. El campo real es 0.37 T a 0.60 m, 0.21 T a 0.65 m y 0.04 T a 0.75 m, y desvía los protones ≈50 mrad. El mapa principal no cambia.
- **Resultados (nivel ×1.0103, puerto 1):**
  - región central: 252/800 (31.5 %), frente al 35.8 % con el campo sintético de KEY22;
  - **el haz sale de la región central a fase +13…+33°**, y la ventana de la v0.34 está en −30…0°. Lo fija la posición de la fuente (40.3° antes del borde del dee, mapa de IBA); el deslizamiento en la región central es solo +8.6°. Se compensa con la corriente del imán: el óptimo pasa de ×1.0098 a **×1.0103** (`coil.escala_iso`; NI ≈ 119.9 kA·v);
  - de la salida de la región central a la lámina llega el 55 % (140/252); con el 15 % de pérdida por gas, el 47 %. TUD03 da un 81 % medido;
  - en la lámina: T = 17.91 ± 0.22 MeV; la lámina se reajusta +1.2° para centrar el haz en el blanco;
  - **en el blanco: σx 7.2 mm, σz 10.5 mm**; dentro de la ventana Ø20 llega el 47 % y por el colimador Ø9 pasa el 20 %;
  - **balance respecto a la corriente DC de la fuente:** 10.5 % sale de la región central, 5.0 % llega a la lámina con gas, 2.3 % entra en la ventana del blanco y 1.0 % pasa el colimador.
- **Diagnóstico de la pérdida principal:** 75 de las 97 pérdidas verticales ocurren entre 14 y 18 MeV. Ahí **νr − 2νz cruza cero dos veces (resonancia de Walkinshaw)**, porque el shim de v0.34 bajó νz a 0.53 a 15–16 MeV. Con el campo medido de KEY22, νz = 0.58 y νr − 2νz se queda en −0.07: no cruza.
  - **Siguiente paso:** repetir el ajuste de v0.34 con una cota νz ≥ 0.57 en 14–16.5 MeV. El jacobiano de νz ya existe en los casos `j2_*`.
  - La mancha en el blanco, mayor que la real, sale del mismo exceso vertical, más el borde del campo (σz' de 6 mrad en la lámina a 13 mrad a 0.6 m).
- `validation/fuente_blanco.{png,json}`.
- **Paquete de ingesta para Radiopharma-OS** (`python herramientas/paquete_rpos.py` → `out/rpos_v035/`):
  - manifiesto de 565 artefactos (STEP, h5m, GDML, STL, 3MF/GLB, Molflow, mapas npz, curvas JSON, parámetros, telemetría), cada uno con md5 y metadatos (formato, unidades, generador, versión);
  - verificación cruzada contra los supuestos KIUBE-V04, el inventario F75, el SRG y KIUBE-OM (22 filas);
  - registro de re-etiquetado append-only e impacto en cadena;
  - telemetría simulada declarada (`origen = simulado_declarado`).
  - Gap efectivo de la cámara: `herramientas/gap_efectivo.py` → `validation/gap_efectivo_camara.json`: 102.6 mm entre caras de polo (V04: 132).

## Qué cambió en v0.34 (ajuste fino de isocronía, adoptado en el CAD)

- **Ajuste fino** (`fem/shim_iso2.py`): 10 perillas Δα por lado entre 100 y 500 mm (`pole.shim_r`, más densas en 400–500 mm) y el borde del escalón 2 (`pole.step2_edge_L`).
  - **Jacobiano medido sobre las órbitas de equilibrio** (11 casos FE perturbados), no sobre ⟨B⟩: el modelo lineal de ⟨B⟩ fallaba en el borde.
  - Mínimos cuadrados sobre el error punto a punto Δf/f(E) **y la fase acumulada del ion central**, con nivel y fase inicial libres, Tikhonov y suavidad del perfil.
  - **Cotas físicas**: el borde de la colina queda a ≥ 5 mm de la caja del dee (Δα ≤ 3.6° a 430–500 mm) y L ≤ 526 mm (el escalón 1 está en 529). La iteración 7, sin cotas, invadía los dees; la 8 es la solución acotada.
- **Resultado (iteración 8, adoptada):** `shim_dalpha = [0, 2.07, 0.69, 0.16, −0.20, 0.53, 0.07, 0.79, 0.03, 3.60, −1.61]`, L = 526 mm.
  - Isocronía de ±0.28 % alrededor del nivel entre 1 y 17.5 MeV. La órbita de 18 MeV cierra (νr 0.81, νz 0.87).
  - Ventana longitudinal 1D del ion central: 96° (KEY22 medido: 104°).
  - En 3D, con betatrón de ±6 mm radial y ±2 mm vertical, llegan a 18 MeV **33 de 70** partículas en el mapa de la iteración 8 (nivel 1.0092; 80 % en las fases −60…−20°). Con el campo medido de KEY22 llegan 47/70.
- **FE oficial del imán** sobre los dominios regenerados: reproduce el mapa de la iteración 8 dentro del ruido de malla (σ 0.07 %, máx. 0.13 %, nivel −0.05 % en el óptimo 3D).
  - En 3D llegan **36/70** con nivel 1.0098 (34/70 con 1.0097 y con 1.0099), frente a 23/70 con 1.0094 y 22/70 con 1.0101: la tolerancia en corriente es de unos ±0.02 % de campo. Por debajo de 1.0092 no llega casi nada.
  - **Nivel de campo de operación: `coil.escala_iso = 1.0098`**: +0.98 % de campo sobre el mapa de NI = 115 kA·v. Con el hierro saturado (d lnB/d lnNI = 0.24) eso exige NI ≈ 119.7 kA·v, unos 20.6 kW de bobinas (IPAC26: 19 kW).
- **CAD:**
  - colinas con huella poligonal (lados a hc ± (hw − asin(o/r) + Δα(r)), paso de 15 mm más los nodos del shim, 80 caras; 0.01 % de volumen respecto al polígono de los FE);
  - **dummy dees**: el corte con la colina ensanchada dejaba astillas sueltas; se conserva solo el sólido principal;
  - **caja del dee `dee.top_z` 76.7 → 74.5 mm**: las colinas más anchas añaden capacidad dee–tierra y bajaban f a 40.10 MHz. Con 74.5, fem_rf da 40.63 MHz (−0.04 %, dentro del rango del sintonizador);
  - el dominio RF usa un polígono de colina de 45 mm (error de forma < 0.05 mm) y `build.py --lod fem_rf` reintenta con malla 0.038 si la factorización se queda sin memoria.
- **Efecto secundario:** en r = 0.45 m el valle baja a 0.42 T en el nivel de operación (S09: 0.45 T, −7 %). El ⟨B⟩ de la banda 0.38–0.47 m coincide con KEY22 en el nivel de operación de cada mapa (+0.03 %). `tests/test_modelo.py::test_calibracion_iman` compara ahora en ese nivel y acepta ±0.04 T en el valle.
- `validation/isocronia.{png,json}` y `validation/aceleracion_3d.png` actualizadas.

## Qué cambió en v0.33 (isocronía del imán y aceleración completa: diagnóstico)

- **Órbitas de equilibrio** en el mapa FE (`fem/orbitas.py`).
  - El mapa se descompone en armónicos 4k con splines en r y se usan las ecuaciones de Gordon con θ como variable independiente.
  - Para cada energía: órbita cerrada, f_rev, νr y νz.
- **Aceleración 3D** (`fem/aceleracion.py`): Boris con la expansión del campo fuera del plano medio y golpes de RF en los 4 bordes de los dees. Incluye el modelo longitudinal del ion central (`deriva_fase`, `mejor_nivel`).
- **Validación con el campo medido de KEY22:** isocronía de ±0.2 % entre 2 y 17 MeV, órbita de 18 MeV a 452 mm y ventana de fases de **104°** (excursión de ±14°, 175 vueltas). Es el comportamiento de la máquina real.
- **Hallazgo: el imán FE de la v0.32 no es isócrono.**
  - Error de frecuencia de −3.7 % a 1 MeV: ⟨B⟩ queda un 3.8 % bajo en r = 0.1 m porque los lados de las colinas, desplazados 10 mm hacia dentro, las estrechan cerca del centro.
  - El campo cae unos 25 mm antes que en KEY22 y no hay órbita cerrada a 18 MeV.
  - El ion central no pasa de 1.3–3.7 MeV con ningún ajuste de la corriente.
- **Ajuste de isocronía (shimming)** (`fem/shim_iso.py`). Perfil Δα(r) de los bordes de colina en 8 radios (`pole.shim_r` y `pole.shim_dalpha`; sin efecto si no se definen) y borde del escalón 2.
  - El interior queda bien ajustado (±0.25 %).
  - **El escalón 2 a 520 mm** (en vez de 512 de T04F2) cierra la órbita de 18 MeV a 449.5 mm (νr 0.89, νz 0.81).
  - Las cuñas del corrector de gradiente son necesarias: al reducirlas, νz cae a 0 cerca de 17 MeV.
  - El mejor caso llega a 18 MeV con una ventana longitudinal de 54° y +1.1 % de nivel, pero en 3D con betatrón solo llega el 2–3 % del haz.
  - **Pendiente:** la zona de 15–18 MeV (r 410–452 mm) necesita perillas más finas que las actuales (cada 70 mm), y una malla FE más fina para bajar el ruido de ±0.1–0.2 %.
  - **No se adopta todavía ningún cambio de polo en el CAD.** `validation/isocronia.{png,json}`.
- **FE del imán más robusto** (`fem/magnet_ngsolve.py`):
  - curva B-H extrapolada hasta H = 1e9 A/m (antes ~8.5 T: el spline dejaba de ser convexo);
  - una búsqueda lineal que trata NaN como «no mejora»;
  - reintento desde una solución lineal (μr = 200) y continuación en la corriente.

## Qué cambió en v0.32 (chimenea y puller en el CAD)

- **Diámetro de la punta de la chimenea** (`validation/region_central_v32.{png,json}`), con el bisel de 5 mm, el centro fijo a 40.3° del borde de la punta del dee y una pared de rendija de 0.4 mm.
  - Más delgada es mejor: Ø5/6/8/10/13 dan 14.7/14.3/13.1/12.1/11.5 % de la corriente DC de la fuente a 6 vueltas.
  - Con Ø13 la ventana de fases aceptadas baja de 90° a 75°, y un puller a 2.8 mm quedaría dentro del cuerpo.
  - **Se mantiene Ø6.** El Ø13 del mapa Opera de IBA se interpreta como el cuerpo de la fuente detrás de la punta.
- **CAD:**
  - **Punta de la chimenea** de Ø6 (|z| ≤ 11 mm), con la rendija de 0.6 × 6 mm hacia la dirección de extracción. Está 40.3° antes del borde de entrada de la punta del dee que alimenta (fuente 0 → dee de 90°, fuente 1 → dee de 270°), igual que en el modelo de la región central.
  - El vástago (Ø10.4) va recto desde la punta hasta el puerto de la fuente.
  - **Mordazas del puller** con el bisel de 5 mm a 45°, a 2.8 mm de la chimenea, fundidas con la punta del dee en `dee_bridge_upper/lower`.
  - Nuevos parámetros en `ion_source`: `chimenea_r`, `ang_antes_borde`, `rendija_wh`, `gap_puller`, `puller_lwh`, `puller_chaflan`.
  - Antes la fuente era un cilindro de Ø13 a 6.7°/186.7° y no había puller.
- **RF:** las mordazas añaden cobre al potencial del dee en el centro. El modo en fase pasa de 40.656 a **40.640 MHz** (−16 kHz), a 6 kHz del objetivo de 40.646 MHz y dentro del rango del sintonizador (±90 kHz); no se recalibra.
- Vista de la región central del CAD en cortes a z = 0, 6 y 10 mm: `validation/cad_region_central.png`.

## Qué cambió en v0.31 (forma de las mordazas del puller)

- **Mordazas parametrizadas** (`fem/central_region.py: mordazas_puller`, compartidas con el modelo del menisco): abertura en dos escalones (`puller_h2`, `puller_d1`) y bisel del borde interior frontal (`puller_chaflan`, `puller_chaflan_ang`).
- **Barrido** (`validation/region_central_v31.{png,json}`):
  - El **bisel a 45°** sube la corriente que llega a 6 vueltas de forma monótona y se satura a 5–6 mm: del 10.3 al **14.3 %** de la corriente DC de la fuente (+38 %).
  - Comparación justa: la transmisión × el ancho de la ventana de fases / 360°, porque el bisel estrecha un poco la ventana aceptada (105 → 90°).
  - El ángulo importa poco (30–60°). Los escalones en la abertura no ayudan.
- **Confirmado con el menisco autoconsistente:** 26.0 → **35.8 ± 2.4 %** en la misma ventana (+37 %).
  - La divergencia vertical a la salida del puller baja a la mitad (219 → 110 mrad rms) y las pérdidas verticales pasan del 40 al 24 %.
  - Lo que queda son sobre todo partículas de fase tardía que la RF devuelve contra la chimenea.
- **Campo máximo** en la región fuente–puller a 40 kV (`fem/central_region_campo_max.py`): 37 → **14 MV/m** (de 4.5 a 1.7 veces Kilpatrick), porque desaparece la arista viva.
- **Adoptado:** `puller_chaflan` = 5 mm a 45° en el preset `kiube`.
- **Pendiente:**
  - El CAD no modela el puller.
  - El cuerpo de la fuente del CAD mide Ø13 mm (mapa Opera de IBA, ±2 mm), mientras que la chimenea del modelo de la región central mide Ø6 mm (EST). Hay que conciliarlos antes de dibujar el puller en el CAD.

## Qué cambió en v0.30 (menisco del plasma de la fuente)

- **Menisco autoconsistente** (`fem/menisco.py`), calculado en un dominio local alrededor de la rendija con una malla de 0.08 mm:
  - el borde del plasma es una frontera móvil; cada punto emite el mínimo entre la j del plasma y la de Child–Langmuir local (φ a 0.1 mm por delante, con la carga espacial de los propios iones);
  - el menisco se mueve hasta que las dos coinciden (en 5 iteraciones, ±10 %);
  - cuasiestático, un menisco por nivel de voltaje: V_dee·cos φ = 0.97 / 0.86 / 0.70 / 0.50;
  - la j del plasma sale de Kuo y Hendry (TRIUMF, Cyclotrons 1992): **144 mA/cm²** en una fuente H⁻ interna de ciclotrón.
- **Forma resultante:**
  - plano en vertical salvo en los extremos de la rendija y **cóncavo a lo ancho** (0.19 mm de flecha);
  - **0.23 mm más hondo** que la cara plana de v0.29 a pleno voltaje y 0.10 mm a medio voltaje;
  - el haz sale del puller muy divergente en vertical (220–265 mrad rms, totalmente correlacionado): esa óptica la fija la lente del puller, no el menisco.
- **Región central con el menisco** (`fem/central_region_menisco.py`): un campo global por nivel, partículas que arrancan en el menisco y malla de 0.08 mm en el canal de la rendija.
  - Transmisión a 6 vueltas: **26.0 ± 2.2 %**, frente a 32.3 ± 2.3 % con el plasma plano.
  - La diferencia es el tiempo de tránsito: con el menisco más hondo, las partículas de fase tardía (+5 a +23°) no cruzan el gap antes de que la RF se invierta y vuelven contra la chimenea.
  - Sensibilidad a la j del plasma: 72 / 144 / 288 mA/cm² dan **22.8 / 26.0 / 29.8 %**.
- **Coherencia:**
  - la rendija de 6 × 0.6 mm emite ≈4.4 mA en las fases útiles;
  - con la ventana de −80 a +40° y el 26 %, llegan ≈0.35–0.4 mA tras 6 vueltas, del orden de lo que el KIUBE necesita para 150–180 µA en el blanco (con la transmisión del 81 % de TUD03T).
- `validation/menisco.{png,json}`.

## Qué cambió en v0.29 (región central: nariz y chimenea)

- **El modelo FE de la región central ahora coincide con el CAD:**
  - el cuello de 3.5 mm que une las puntas (TUD03T);
  - la caja del dee a ±14 mm y 39° desde r = 60 mm, con el escalón punta → caja;
  - el tapón rebajado 15 mm (cara a z = 27; antes tocaba las puntas a z = 12);
  - el límite vertical de pérdida es ±7.8 mm dentro de las puntas y ±11.8 mm fuera.
- **Barridos** (`validation/region_central_v29.{png,json}`), con 400 partículas, 6 vueltas y ±2.3 %:
  - **Giro de la rendija de la chimenea:** la meseta va de 0 a +10° (tangente o algo hacia fuera). Hacia dentro cuesta 6–8 puntos y +20° cuesta 16. Se mantiene **tangente (0°)**.
  - **Radio de la punta** (6–10 mm) y **ángulo de las puntas** (44–56°): sin efecto fuera del error. La forma medida en el mapa de IBA (r 8 mm, 50°) está en la meseta.
  - **Semiabertura de las puntas:** con 7 mm se pierden 6 puntos; con 9–10 mm se ganan +2–3 (≈1σ). Se mantiene 8 mm para no recalibrar la RF.
  - **Cuello entre puntas:** −4.5 ± 3.3 puntos, pero centra mejor la órbita (7.8 frente a 9.2 mm). Se conserva porque es la geometría medida.
- **Diagnóstico** (`fem/central_region_diag_z.py`):
  - El 61 % del haz se pierde por la vertical en el **primer cuarto de vuelta** (mediana 0.4 vueltas, r ≈ 18–25 mm). Son las partículas que salen con |z| > 0.8 mm de la rendija de 6 mm.
  - Es la óptica de extracción de la fuente: el modelo arranca los iones de una superficie plana, sin el menisco del plasma, que en la realidad enfoca. La transmisión absoluta (≈31 %) es una **cota inferior**; las comparaciones relativas siguen valiendo.
  - La geometría de la nariz y de la chimenea ya no es el factor limitante.

## Qué cambió en v0.28 (techo del búnker, preliminar)

- **Cálculo del techo del búnker del proyecto** con un cuarto de aire acondicionado encima (ocupación mínima): `mc_benchmark/bunker_techo` (Geant4) y `validation/techo_bunker.py`.
  - Búnker de 5.9 × 4.7 × 2.5 m con muros de 2.11 m, ciclotrón centrado, puertas de HDPE, hormigón de 2.35 g/cm³.
  - Reducción de varianza hacia el techo; curva dosis–espesor en una sola corrida, comprobada con techos reales de 1.5 y 2.0 m.
- **Criterio (México, RGSR):** 1 mSv/año (Art. 38) con T = 1/40, y menos de 50 µSv/h durante la irradiación (Art. 113).
- **Resultado:** TVL de 24.6 cm.
  - Para 50–150 k µAh/año hacen falta **1.55–1.74 m** escalando al SRG (1.27–1.46 m con el modelo directo).
  - El mínimo del SRG (1.70 m) cumple, salvo con 150 k µAh/año y T = 1/20.
  - Pendiente de los cortes: altura real, espesor, posición del ciclotrón, penetraciones de ductos y tapa de izado.

## Qué cambió en v0.27 (blancos en la cara del yugo y puertas blindadas)

- **Blancos en la cara del yugo:** en VE24 Fig. 3 (ejes en metros) los blancos están fuera de la cara del yugo (±0.88–0.98 m), y la foto de VUB 2024 muestra las puertas blindadas en esa cara. Los 8 blancos pasan de r = 590 mm (EST, junto a la cámara) a la cara + 10 mm, con un tubo de haz por la escotadura del puerto. Los puertos #1/#2 de VE24 caen a ~25 mm de los del modelo. **La fuente de neutrones pasa de r = 597 a 969 mm** (~37 cm hacia fuera), lo que importa para el blindaje.
- **Puertas blindadas Nirta** de **HDPE de 18 cm** (ficha IBA citada por VUB) en los blancos F18 del proyecto (puertos 1 y 5). Nuevo material `hdpe`.
- **Ventana del blanco F18: Havar de 35 µm** (IBA, Cyclotrons2013), antes 50.
- **Validación de activación** de un KIUBE real (VUB-UZ Brussel): Eu-152 y Co-60 en el hormigón hasta 1 m frente a los blancos, y el cálculo de IBA con y sin puertas → `validation/activacion_vub_kiube.json`.
- La envolvente de maniobra del SRG se comprueba sin los blancos, que se montan en sitio.
- **Benchmark Geant4 frente al SRG repetido** con la nueva fuente (sin puertas, como las Tablas 8–9): dentro del búnker el modelo pasa de 0.23–0.54 a **0.51–1.22 del SRG**; tras 2–2.3 m de hormigón, de 0.007–0.02 a 0.06–0.08; TVL de 28.6 cm frente a 33.9 cm (ajuste ponderado por 1/σ en todos los espesores). `mc_benchmark/correr_caso.sh` corre un caso suelto (espesor, semilla).
- **Puertas de HDPE en Geant4:** la fluencia en el muro frontal baja ÷1.85–2.7 (IBA: activación ÷2.3) y el perfil en profundidad coincide con las varillas medidas en el KIUBE de Bruselas → `validation/benchmark_puertas.{json,png}`.
- **STL HD:** se descartan las láminas de espesor cero que dejaban los tubos de haz partidos por el plano medio; el modelo completo vuelve a ser un solo sólido (18 puntales de impresión).

## Qué cambió en v0.26 (revisión de cotas estimadas)

- **Informe `docs/COTAS_EST.md`:** las 36 cotas que siguen como EST, ordenadas por impacto, con qué afectan y qué documento las cerraría.
- **`dee.angle` 38 → 39°**, medido en el render en planta de TUD03T (diap. 9) rectificado con los insertos de polo (39.1 ± 0.6°). RF recalibrada: caja del dee a 76.7 mm, **40.656 MHz**.
- **Reclasificadas** 4 cotas del vacío que no son estimaciones: constantes físicas (`CONST`) y opciones del modelo (`MODELO`).

## Qué cambió en v0.25 (rendija de la fuente)

- **Rendija de la chimenea** barrida de 3 a 8 mm: la transmisión baja de 55 a 25 %, pero la corriente útil (∝ altura × transmisión, porque la PIG emite en proporción al área) es máxima con **6–8 mm**. La rendija se queda en 6 mm (ahora CALC). Figura en `validation/region_central_rendija.png`.
- La pérdida vertical que queda (≈ 62 %) ocurre en la primera vuelta y viene de la altura de salida de los iones: es el compromiso entre emitancia y corriente de la fuente.

## Qué cambió en v0.24 (forma del puller y abertura central)

- **Puller optimizado** (región central, 400 partículas, 6 vueltas): **abertura vertical de 10 mm** (antes 5, EST), ancho de 7 mm y profundidad de 3 mm. La transmisión del haz sube de 19.5 a **≈ 35 %**: con 5 mm el 37 % del haz chocaba con las mordazas del puller nada más salir de la fuente. Figura en `validation/region_central_puller.png`.
- **Puntas de los dees a ±8 mm del plano medio** (`dee.punta_z0`; antes ±14 como la caja). Con ±14 la región central no acelera: aceptancia de 20–30° y 82–90 % del haz perdido en la primera vuelta. Un escalón une cada punta con su caja.
- **RF recalibrada:** caja del dee a 78.7 mm → 40.652 MHz (+0.01 %).
- **Corrección de v0.23:** el techo de transmisión no es el 29 % (105°/360°) sino ≈ 87 %, porque las partículas se lanzan dentro de la ventana de fase aceptada ±10°.

## Qué cambió en v0.23 (gap fuente–puller)

- **Barrido del gap chimenea–puller** en la región central (malla fina, 6 vueltas, 400 partículas, fuente en la posición de IBA): 2.0 / 2.3 / 2.6 / 2.8 / 3.2 mm → 14.8 / 18.2 / 21.5 / 19.5 / 19.2 % de transmisión (±2 %).
  - El óptimo está en **2.8 mm**: con menos gap el haz choca contra el puller, con más se pierde en vertical.
  - Con 2.8 mm la fuente queda a 40.3° del borde del dee, como en el mapa de IBA.
  - El campo en el gap es de 14.3 MV/m (1.7 Kp).
  - El gap del C18/9 escalado a 40 kV (2.0 mm) transmite un 7 % menos.
- Preset `kiube` de `fem/central_region.py`: `gap_sp` = 2.8 (CALC). Nuevo `fem/central_region_barrido_gap.py` para barridos rápidos. Figura en `validation/region_central_gap.png`.

## Qué cambió en v0.22 (región central del mapa Opera de IBA)

- **Puntas de los dees y cuello** medidos en el mapa Opera de la región central (TUD03T, diap. 5; la figura tiene escalas distintas en x, 10.3 px/mm, y en y, 7.8 px/mm):
  - puntas de 50° (bordes radiales a ±25°) desde r = 8 mm hasta la caja del dee;
  - unidas por un cuello de 3.5 mm, desplazado 1.25 mm del eje;
  - placas de 4 mm (EST): el campo hasta el tapón baja a 4.4 MV/m (0.5 Kp).

  Sustituye a la placa de 36 mm de v0.19 (`dee.punta_*`, `dee.puente_*`). Vista en `validation/region_central_v22.png`.
- **Fuentes gemelas a r = 12.5 mm** (antes 18, EST), con un cuerpo de chimenea de Ø13. En el preset `kiube` de `fem/central_region.py`, la chimenea queda 39.4° antes del borde de entrada de su dee; el mapa de IBA da 40°.
- **RF recalibrada:** caja del dee a 77.9 mm → **40.668 MHz** (+0.06 %). La contrafase queda en 110 MHz y el acoplador sigue con un lazo de 30 mm para 22.2 kW.

## Qué cambió en v0.21 (imán cerrado con 19 kW y las patentes de IBA)

- **Patentes de IBA de 2016** (`PAT16`: US9907153, US10064264, US10278277, más US9961757 ya registrada). Dan rangos, no cotas; `validation/patentes_iba.py` comprueba 20 reglas.
- **Referencia de potencia: 19 kW de bobinas** (IPAC26, decisión del usuario) → `NI_total` = 115 kA·v.
- **Imán recalibrado dentro de las patentes:**
  - valle a 100 mm (Gv = 200; en v0.18 era 95, fuera de la reivindicación 19);
  - **corrector de gradiente en cuña de 90°**, de 78 × 39 mm (`pole.gc_forma`);
  - tapón a −15 mm.
- **Resultado a 115 kA·v = 19 kW:** ⟨B⟩ = 1.3253 T (IBA 1.3253), flutter 0.389 (0.394), joroba central 1.030 (1.033), colina/valle 2.22/0.426 T. El acero no se toca (J_s 2.12 T) y la malla fina cambia ⟨B⟩ un −0.05 %.
- **RF recalibrada:** caja del dee a 75.6 mm → 40.626 MHz (−0.05 %).

## Qué cambió en v0.20 (6 documentos de RF y la presentación del KIUBE)

- **Presentación TUD03 del KIUBE** (`TUD03T`):
  - **El mapa Opera de la región central confirma que los dees están unidos en el centro**, por un cuello estrecho, como suponía el puente de v0.19. También muestra puntas de dee a r ≈ 10 mm y chimeneas a r ≈ 12 mm; quedan anotadas, sin aplicar todavía.
  - **Vacío calibrado** con el bombeo medido (4 × VHS-6): el acero dulce del yugo desgasifica 30 veces más que el metal limpio (`outgassing_factor.pared_acero` = 30). Antes el modelo era 15–40× optimista. El vacío base pasa a 2.7e-7 mbar. Prueba nueva: `test_bombeo_frente_a_tud03`.
  - **Discrepancia del imán**: «< 11 kW» frente a los 19.4 kW que pide el FE calibrado (ver `docs/FUENTES.md`).
- **KIRAMS-30** (`KIR30`): el modelo de línea de transmisión valida el FEM de RF (C_dee 154 frente a 167 pF, 8 %).
- **Pick-up RFT-30, curso CAS de Brandenburg, C70 y C400 de IBA:** referencias. Confirman el criterio de sobreacoplamiento sin haz, el comportamiento inductivo del pistón y que IBA usa tetrodos y CST.

## Qué cambió en v0.19 (acoplamiento RF, acoplador y sintonizador)

- **Dos papers revisados** (`docs/FUENTES.md`):
  - *Instruments* 8 (2024) 38 es **VE24**, ya incorporado. Aporta las cuñas definidas (13/14/15/18 MeV): nuevo `extraction.ve_energias_cuna`.
  - Cyclone 11 (`CY11`, TOSCA): su sensibilidad d ln⟨B⟩/d ln NI de 0.21 confirma la saturación del FE calibrado (0.24). También da el carrusel de diseño C18/9: `fork_len` 25 mm y 2 láminas.
- **Puente central entre dees** (`dee.puente_w` = 36, EST). Un solo amplificador y un solo cable de 50 Ω (SRG) exigen dees acoplados. El puente lleva la contrafase de 40.7 a 125 MHz y deja el modo acelerador aislado y equilibrado (V1/V0 = 1.00).
- **Sintonizador** (pistón de Ø40 en la carcasa del dee 1): rango de ±90 kHz por pistón alrededor de la posición nominal, unas 4 veces el ancho de banda de la cavidad cargada.
- **Acoplador** (lazo en la carcasa del dee 0, `fem/rf_acoplos.py`): **30 × 70 mm** para acoplamiento crítico con **21.8 kW** (18.5 kW de paredes con Q real = 0.65 × cobre + 3.24 kW de haz a 180 µA). Q cargado 1683 y ancho de banda de 24 kHz. Coherente con los 12 kW de calor del amplificador en el SRG (eficiencia del 64 %).
- **RF recalibrada:** caja del dee a 72.6 mm → **40.683 MHz** (+0.09 %), Q de cobre 5658, 12.0 kW a 40 kV. La malla de RF pasa a maxh 0.035 m con gradación 0.55: la anterior (0.045 y 0.4) subestimaba f un 0.3 % y ya no cabía en memoria con el puente y el pistón.
- **Clasificación de modos robusta** en `fem/rf_kiube.py`: en fase, contrafase u otro, según V1/V0. Figura en `validation/rf_acoplamiento.png`. Detalle en `fem/README_RF.md`.

## Qué cambió en v0.18 (calibración del imán)

- **Diagnóstico nuevo:** `fem/magnet_ngsolve.py --diag` reparte ∮H·dl = NI/2 por tramos del circuito magnético. Con la curva 1010 genérica el hierro se llevaba el 75 % de la fuerza magnetomotriz, sobre todo en las puntas de polo (2.35 T) y la tapa A (2.22 T). Esa curva satura en J_s ≈ 2.0 T, y el acero de bajo C de los ciclotrones llega a 2.10–2.15 T.
- **Recalibración** con KEY22 y S09 (`validation/calibrar_iman.py`, tabla y figura `validation/calibracion_iman.png`):

  | Parámetro | Antes | Ahora | Por qué |
  |---|---|---|---|
  | Curva B-H del yugo | 1010 genérica (J_s 2.0) | `bh_acero_bajoC_kiube` (J_s 2.12 T) | acero de bajo C |
  | Tapa A | 166 mm | 181 mm | extremo superior de T04F1 (±15) |
  | Suelo del valle | 110 mm | 95 mm | flutter de KEY22 |
  | Rebaje del tapón | 17 mm | 16 mm | joroba central de KEY22 |
  | Bobina, z inferior | 74 mm | 59 mm | sigue a la tapa |
  | NI que pide el FE | 170 kA·v | **116 kA·v** | potencia: 105 (100–115) |

- **Resultado a 116 kA·v:** ⟨B⟩ = 1.3240 T (IBA 1.3253), flutter 0.381 (0.394), colina/valle 2.20/0.446 T (S09: 2.15/0.45), joroba central 1.026 (1.033). La malla de 6 mm cambia ⟨B⟩ un −0.07 %.
- **RF recalibrada** por el nuevo valle: la caja del dee baja de 86 a 75 mm (20 mm bajo el suelo del valle, 2.0 MV/m a 40 kV) y la frecuencia queda en 40.70 MHz (+0.14 %). Densidad de `copper_dee_caja`: 3.174.
- **Geometría:** los carruseles de strippers bajan para no chocar con la bobina. La masa total pasa a 17.8 t (SRG: 18 t).
- **Opera-3D/TOSCA:** el paquete lleva la curva calibrada y NI = 116 kA·v.
- **Corrección de v0.17:** la "eficiencia BATAN" de 0.84 no se puede trasladar al valle profundo del KIUBE (sale 0.58), así que los "≈100 kA·v" deducidos de ella no valían. Queda como cota en `docs/FUENTES.md`.

## Qué cambió en v0.17

- **Uso personal:** se retiran las advertencias de confidencialidad (los documentos son del usuario, con derechos de uso).
- **Tres papers revisados** (`docs/FUENTES.md`):
  - `BATAN13`: queda como chequeo del circuito magnético. Su eficiencia es 0.84 frente a 0.50 en nuestro FE, lo que apunta a unos 100 kA·v y no a 170 k.
  - `AUT15`: solo como referencia del acero AISI 1010.
  - `KEY22`: confirma la topología RF.
- **Paquete Opera-3D/TOSCA regenerado** en `out/opera3d` (`python fem/opera_export.py`), con NI 105 kA·v.

## Qué cambió en v0.16

- **Tres documentos revisados** (`docs/FUENTES.md`).
  - **Se incorpora la clase JUAS 2022 de IBA** (`JUAS22`, pública).
  - Se desechan el S2C2 (sincrociclotrón de protonterapia) y la fuente ECR de 28 GHz: no aportan nada al KIUBE.
- **Dummy dees a tierra** (fotos del diseño de valle profundo, JUAS22 p.44–45): placas de cobre apoyadas en los bordes de las colinas, a 6 mm de los bordes del dee (EST). Forman los 4 gaps aceleradores. Ahora hay 8 piezas más, sin solapes.
- **RF recalibrada con los dummy dees.**
  - La caja del dee baja de 97.6 a 86 mm, a 24 mm del suelo del valle (1.7 MV/m a 40 kV, más margen de tensión).
  - Resultado: 40.70 MHz, Q (Cu) de 5895 y 11.7 kW a 40 kV.
  - La alternativa de dee plano que sugiere la foto (placas con dummy dees a 3 mm) da 46 MHz y queda registrada.
- **Región central del C18/9:** gap chimenea–puller de 2.6 → 1.5 mm (JUAS22). El KIUBE mantiene 2.6 mm (EST).
- **Datos de instalación con respaldo público** (JUAS22: plano medio a 1.12 m, altura de 1.85 m, 1.95 × 1.95 m, 17 t). La versión pública ya usa el plano medio real (1550) en vez del aproximado (1600).
- Todo regenerado:
  - 10 LOD sin solapes; STL HD con 62 cuerpos estancos;
  - pruebas: 16/16;
  - SRG: 9/9, con masa de 17 523 kg (SRG 18 t, JUAS22 17 t).

## Qué cambió en v0.15 (consolidación)

- **Todas las salidas regeneradas** con la geometría actual: neutronics, fem_struct, print, fem_mag, vacuum, render, print3d, mc, stl_hd y fem_rf. Ninguna tiene solapes ni sólidos inválidos.
  - `stl_hd`: 62 cuerpos estancos y manifold, sin autointersecciones; el completo es un solo sólido.
  - `mc`: 100 volúmenes, con relación de volumen del teselado entre 0.994 y 1.001.
- **La caja del dee es maciza con cobre homogeneizado** (`copper_dee_caja`, ρ = 2.574). La caja hueca cerrada rompía las mallas de vacío y Monte Carlo. La superficie de RF es la misma, y la masa también (17 521 kg).
- **Batería de pruebas automáticas** (`tests/`, pytest): 16 pruebas, todas en verde.
  - Rápidas: fuentes declaradas, materiales, validez y solapes en 3 LOD, masa, cotas de instalación, TPMC contra Clausing, paleta, STL HD, 3MF estricto (lib3mf) y validación SRG.
  - Lentas: RF contra soluciones exactas y frecuencia del modo acelerador.
  - **Ya encontró un fallo real:** el 3MF de color tenía 2 triángulos degenerados (en el colector de latón). Ahora se filtran.
- **Versión pública** (`herramientas/version_publica.py`).
  - Genera `dist/kiube_publico/`: las cotas SRG/PROP vuelven a sus valores públicos y desaparecen los bloques `proyecto`, `sitio` y `fuente_rad`.
  - Excluye el benchmark de blindaje y limpia README y comentarios.
  - Después audita que no queden referencias (0) y corre las pruebas en la copia (11 pasan; 3 se omiten por falta de salidas o del SRG).
- No se pudieron repetir aquí, porque el reinicio del entorno borró los binarios:
  - la comparación con Molflow (falta `molflowCLI`; el motor TPMC propio sí pasa las pruebas contra Clausing);
  - los renders de Blender.

## Qué cambió en v0.14

- **Modelo de RF** (`python build.py --lod fem_rf`; detalle en `fem/README_RF.md`).
  - Modos propios del resonador con NGSolve (Nédélec de orden 2). Verificado contra soluciones exactas: cavidad TM010, su Q y coaxial λ/4, con errores de 10⁻⁶.
  - **Con la geometría estimada anterior el resonador daba 79 MHz, no 40.646.** El dee de 2 placas finas no tenía capacidad suficiente.
  - Ahora el dee es una caja hueca de cobre cuya cara superior está a 97.6 mm, 12.4 mm bajo el suelo del valle. El stem pasa de Ø90 a Ø50.
  - Resultado calibrado: **40.664 MHz** (+0.04 %), Q (cobre) de 5870 y 12 kW de pérdidas a 40 kV.
  - El par de cotas calibrado no es único; el mapa f(altura del dee, radio del stem) está en `out/fem_rf/barrido.json`.
- La masa del modelo queda en 17 521 kg (SRG: 18 000).
- Las salidas de vacío, Monte Carlo e impresión deben regenerarse con los dees nuevos.

## Qué cambió en v0.13

- **Benchmark de blindaje con Geant4 frente al SRG** (`mc_benchmark/`, `validation/benchmark_srg.py`; detalle en `mc_benchmark/README_BENCHMARK_SRG.md`).
  - Montaje: el ciclotrón completo en el búnker estándar del SRG, doble haz F18 a 2 × 100 µA y H*(10) según ICRP 74.
  - Reducción de varianza por importancia geométrica en el muro.
  - Dentro del búnker, el modelo da 0.23–0.54 veces el SRG. A 1 m del blanco coincide con el cálculo analítico isótropo; la diferencia apunta a la emisión hacia adelante que el modelo no tiene.
  - **Tras 2.0–2.3 m de hormigón da 50–150 veces menos que el SRG** (TVL de 25.6 cm frente a 33.5 cm).
  - Con la mitad de hidrógeno en el hormigón la dosis sube unas 10 veces: el agua del hormigón es la incertidumbre dominante.
  - **Para el diseño mandan el SRG y el Monte Carlo de IBA.**
- Geant4 11.4.3 con G4NDL 4.7.1 desde conda-forge (el instalador por pip no traía los datos).

## Qué cambió en v0.12

- **Tres documentos nuevos** (registro, confidencialidad y discrepancias en `docs/FUENTES.md`):
  - el Site Requirement Guide de IBA (Rev. G, 2023);
  - el folleto de 2017;
  - la cotización del proyecto Zumaya (2026).
- **Cotas de instalación (SRG).**
  - Plano medio a 1120 mm del piso (`floor_to_top` pasa de 1600 a 1550).
  - Cáncamo central M48 con altura total de 1850 mm.
  - Pies a ±525 mm (pilar de 220 y pie de 350).
  - Actuador de la fuente acotado a la envolvente de maniobra (2050 mm).
- **Configuración del proyecto (PROP).**
  - KIUBE 180 con doble haz y 2 fuentes.
  - Blancos: 2 F18-700 (3.6 ml) y 1 13N (2.9 ml), en los puertos 1, 5 y 3 (asignación EST).
  - 4 strippers de serie.
  - La profundidad de cada cavidad sale ahora de su volumen (`target.volumen_ml`).
- **Bloques nuevos en `params.yaml`:**
  - `proyecto`;
  - `sitio` (agua, potencia, gases, cable RF, búnker, armarios);
  - `fuente_rad` (n/p y γ/p del SRG y dosis de referencia para validar el Monte Carlo).
- **Monte Carlo.**
  - El término fuente toma 180 µA y doble haz del proyecto.
  - Normalización `BON21` (por defecto) o `SRG`.
  - Hormigón de 2.35 g/cm³ (mínimo del SRG).
  - Techo del búnker de 1700 mm (antes 1500, por debajo del rango del SRG).
- **Autoblindaje.** Envolvente SRG con pared de 410 mm (calculada para que cada mitad pese 48.1 t), en el material homogeneizado `selfshield_magnetita` (acero, magnetita y agua).
- **Validación** (`validation/validar_srg.py`): 9 de 9 comprobaciones bien. **La masa del modelo es 17 575 kg frente a 18 000 del SRG (−2.4 %)**, una comprobación independiente del CAD.

## Qué cambió en v0.11

- **Servicios de la máquina real** (`kiube/auxiliares.py`, bloque `auxiliares` de `params.yaml`). Son 43 piezas cosméticas (solo LOD print y render), con cotas medidas en la foto VE24 p.7 tomando la altura del yugo como escala (etiqueta `F7`, ±20 %):
  - colector de latón DN40 con 8 ramales y válvulas de bola (manilla roja);
  - mangueras azules a las bobinas y 4 transmisores azules;
  - acometida negra de EPDM;
  - caja eléctrica;
  - manguera neumática en espiral;
  - cadena portacables sobre la tapa.
- Se activan o desactivan con `auxiliares.activo` (la máquina «limpia» del render oficial) y se cambian de cara con `cara_servicio_az`.
- Sin solapes con el resto del modelo, todo estanco y manifold. Neutrónica, FEM y vacío no cambian (piezas solo cosméticas).
- Nuevos cuerpos imprimibles: `servicios_agua` (6 colores, partido en `por_color/`), `cadena_portacables`, `caja_electrica` y `manguera_neumatica`. Ahora son 65 cuerpos.
- Materiales nuevos en `materials.yaml` (latón, PU, EPDM, armario y transmisor homogeneizados) y colores nuevos en la paleta (manguera, válvula, transmisor, caja).
- Render: vista `--vista servicios`, con el encuadre de la foto; el estudio gira con la cámara.
- `stl_hd`: ángulo mínimo por cuerpo para tubos finos. La flecha sigue igual; la espiral bajó de 1.6 M a 75 k triángulos.

## Qué cambió en v0.10

- **Paleta de colores medida en las fichas técnicas** (`palette.yaml`, `kiube/palette.py`). Es una sola fuente para el render (Blender), el 3MF/GLB/STL de impresión y los visores.
  - Muestreo reproducible de regiones definidas de VE24 p.1 (render oficial), VE24 p.7 (foto de la máquina) y KEY22 p.4 (foto del interior).
  - Las imágenes CMYK pasan a sRGB con el perfil FOGRA39 (sin ese paso el cobre salía violeta).
  - Corrección de blanco y negro en luz lineal para quitar el viraje azul de las fotos.
  - Metales con la reflectancia F0 de tablas PBR.
  - Cada color lleva fuente, RAL aproximado (ΔE76), filamento sugerido y nivel de confianza (alta / media / baja = estimado).
- Asignación pieza → color común (`palette.VIS`). La foto real corrige el render en tres puntos: chaflanes del yugo gris azulado, gatos de elevación pintados gris azulado y marcos de puerto gris claro. Además, el colector de agua es de latón.
- Los insertos dejan el azul didáctico de v0.9: son acero mecanizado (sigue disponible con `--paleta didactico`).
- `stl_hd`:
  - color por cara (cada cara toma el material de la pieza CAD de origen; los chaflanes del yugo, la banda gris);
  - 3MF con color propio (núcleo 3MF, *basematerials*; los cuerpos multicolor van como objeto con partes);
  - GLB con materiales PBR, `por_color/` y `lista_colores.csv`.
  - **Corrección**: el 3MF de v0.9 no llevaba color (trimesh no lo exporta), aunque el README lo decía.

## Qué cambió en v0.9

- **STL de alta definición para impresión 3D** (`python build.py --lod stl_hd`, ver `blender/README_STL_HD.md`).
  - Los 61 cuerpos por pieza y el completo como un solo sólido (unión booleana con manifold3d, con puentes de impresión declarados). También van la versión multicuerpo y un 3MF con colores.
  - Teselado desde el CAD exacto con flecha absoluta. Antes CadQuery usaba tolerancia relativa: la "flecha en mm" anterior no era real.
  - Validación por cuerpo: estanco, variedad (manifold3d) y auto-intersecciones con verificación exacta. Remallado isotrópico opcional con MeshLab (pymeshlab), con su desviación medida.
  - Variantes 1:10 (alta / ultra), 1:5 y 1:1. Compatible con `blender_prep_impresion.py` (pines, placas, cuadrantes).
- Blancos: un solo taladro cuando la ventana y la cavidad tienen la misma Ø. Elimina una cara de ancho cero que plegaba la malla.

## Qué cambió en v0.8

- **Molflow+ 2.11.1 integrado** (`kiube/molflow.py`).
  - Escribe el XML de Molflow directamente; ya no hay que importar el STL en la GUI. Ejecuta molflowCLI y lee los resultados.
  - `python build.py --lod vacuum --molflow auto` corre Molflow y usa su perfil de presión para las pérdidas de H⁻.
  - Comparado con el TPMC propio (`validation/molflow_vs_tpmc.py`):
    - Clausing: los tres coinciden dentro de 1σ.
    - Pump-down del tutorial: 4.943 frente a 4.936·10⁻⁵ mbar.
    - Cámara del KIUBE: +0.5 % (H₂O) y −1.0 % (H₂).
  - Unidades del XML verificadas: desgasificación en Pa·m³/s por faceta, `desType = 2` = coseno.
- **Región central** (Zaremba & Kleeven, CAS 2017, arXiv:1804.08961).
  - Métrica de ajuste longitudinal (puntos V = 0 sobre ejes cada 180°/h) y guía de fase −40…−10°.
  - Azimut de las fuentes optimizado (+4°): desviación de 12.9° a 7.6° para H⁻.
  - Recorte paramétrico de la chimenea D⁻ para el paso del H⁻.
  - Figura de fuerza vertical eléctrica (equivalente a la Fig. 40).
  - Chimeneas de Cu-W en el render.
- **Imán, Ec. 31 de CAS:** con la caída en el hierro del C18/9 publicado (60 %), el KIUBE pide ≈ 101 kA·v, coherente con los 105 de la potencia. La curva B-H genérica (170) es la que sobra.
- **Blanco líquido de ⁶⁸Ga** (tipo `Ga68`, VE24).
  - ⁶⁸Zn(NO₃)₂ con 33 mg/mL de ⁶⁸Zn y ventana de Nb de 75 µm a 13 MeV o de 250 µm a 18 MeV.
  - El material admite enriquecimiento isotópico en OpenMC.
- **Pendiente con motor comercial:** Opera-3D (paquete listo; NGSolve lo sustituye). SynRad+/SynradG4: cerrado, no aplica.

## Qué cambió en v0.7

- **Imán en FE 3D, equivalente a Opera-3D** (`fem/`, ver `fem/README_FEM.md`).
  - NGSolve no lineal con curva B-H, sector de 90° periódico.
  - Paquete para Opera-3D con geometría, materiales, conductor, simetrías, malla y exportación, más un comparador.
  - Contra IBA (KEY22 Fig. 3) se fijaron dos cotas: **valle a 110 mm** del plano medio (antes 66) y **tapón F rebajado 17 mm**.
  - Ahora coinciden el perfil de flutter y la joroba central: 1.04 frente a 1.03.
  - Hallazgo: con la curva 1010 representativa hacen falta ≈ 170 kA·v, contra ≈ 105 kA·v que permiten 19 kW. La curva B-H real es la siguiente pieza.
- **Región central 3D del Cyclone 18/9** (`fem/central_region.py --preset c18_9`).
  - Laplace 3D con malla graduada: 0.2 mm en el gap fuente–puller, hasta 18 mm en la estructura completa de los dees (r = 480 mm); 287k tetraedros, 422k gdl.
  - Fuentes H⁻ y D⁻ con puller, dees de 30° y dummy dees; todo parametrizado y con optimizador de la posición de las fuentes.
  - Resultado del optimizador: H⁻ a 12.8 mm y D⁻ a 24 mm, la D⁻ más afuera como en CAS 2017 Fig. 37.
  - Órbitas con física real: Lorentz relativista con RF, puntos V = 0 y V máx, aceptancia en fase de 80° (H⁻) y 65° (D⁻).
  - Render 3D en Blender. Paquete Opera-3D y lectura de la tabla de E que exportes de Opera (`--campoE`).
- **Vacío** (`out/vacuum/README_MOLFLOW.md`).
  - STL de la cámara con desgasificación por material del CERN: metales 3·10⁻⁹/t[h] (Chiggiato, CAS 2017) y Viton 1·10⁻⁷ a 10 h (Chiggiato, CAS 2012).
  - TPMC propio (misma física que Molflow) validado con el pump-down del tutorial de Molflow, en modo dependiente del tiempo: 4.94·10⁻⁵ frente a 5·10⁻⁵ mbar, τ 1.56 frente a 1.57 ms. También con Clausing (±0.5 %) y el caso CERN de 1 m².
- **SynRad+ / SynradG4:** no se implementan (`validation/synrad_evaluacion.json`). Un H⁻ de 18 MeV radia ~1·10⁻¹⁰ eV por vuelta, con fotones de ~1·10⁻⁶ eV, 9 órdenes por debajo del corte de Geant4.
- **Render en GPU:** `blender/render_8k_gpu.bat` hace la exterior y la abierta en 8K y la región central en 4K.
- Masa total: 17.6 t (ficha: 18 t); todos los LOD reconstruidos sin solapes (Geant4: 0).

## Qué cambió en v0.6

Nuevos documentos: do Carmo et al., *Instruments* 8 (2024) 38 (energía variable, **VE24**), Nuttens et al., Cyclotrons2022 WEAO03 (Cyclone KEY, **KEY22**) y Herrod et al., IPAC2026 MOP7178 (**IPAC26**).

- **Fuentes de iones radiales.** VE24 Fig. 1 y el KEY indican inserción radial por el plano medio, y el render oficial no tiene caja arriba.
  - Desaparecen las perforaciones verticales de A y F y los posicionadores superior e inferior.
  - Cada fuente es ahora una chimenea en el gap del tapón, un vástago delgado y un tubo portador. El tubo cruza el valle, la cámara con un pasamuros sin fuga, el espacio entre bobinas y la abertura (2) del yugo.
  - Ejes a 6.7° y 186.7°, sacados de la posición del dispositivo «#IS» en VE24 Fig. 3.
  - Por fuera llevan la pieza blanca vertical y el actuador horizontal del render.
  - Las aberturas (2) de T04F2 eran la entrada de las fuentes, no conexiones de bobina.
- **Soportes de blanco eliminados.** La pieza blanca del render pertenece a la fuente; las ventanas de los puertos quedan desnudas.
- **Energía variable 13–18 MeV** (`extraction.ve_energy_MeV`, `ve_ports`). Se calculan por trazado la horquilla larga del stripper y la brida en cuña, y se modela la brida.
  - A 13 MeV la horquilla entra 70.6 mm más y el azimut del stripper cambia solo −1°.
  - La cuña resulta de −21°, −18°, −13°, −8.6° y −4.3° para 13, 14, 15, 16 y 17 MeV.
  - Solo cabe en los puertos del sistema 2 (desplazamiento de 107 mm frente a 125 mm de semiancho de la abertura), igual que VE24: solo dos puertos, #4 y #8.
- **Bobinas (IPAC26).** El cobre total es 2 t, así que el material homogeneizado queda con 74 % de Cu, 16 % de epoxi-vidrio y 10 % de agua (densidad 7.02).
  - NI = 105 kA·v sale de la potencia de 19 kW: ≈ 16 kW de bobinas y J ≈ 2 A/mm².
  - El ⁶⁰Co de las bobinas viene de ⁶³Cu(n,α).
- **Validación** (`validation/validar_documentos.py`, sin ajustar nada):
  - Trayectorias de 18 y 13 MeV contra las calculadas por IBA (VE24 Fig. 3): diferencia máxima de 2.8°. A r = 820 mm el modelo da 17.2° y 44.5°, frente a 17.7° y 42.9°.
  - ⟨B⟩ hard-edge de 13 316 G frente a 13 253 G de IBA (r = 0.38–0.47 m), y flutter de 0.407 frente a 0.394 (KEY22 Fig. 3, digitalizada en `validation/`).
  - Masa total de 17.87 t frente a 18 t de la ficha.
- **Render en GPU:** `--gpu auto` usa OptiX o CUDA; con 24 GB de VRAM el 8K sale en una pasada.

## Qué cambió en v0.5

- **Exterior completado con la ficha técnica** (Cyclone KIUBE V6, ago-2026, y sus renders):
  - Carcasas RF huecas sobre los valles de los dees, bridas ciegas y 2 bombas difusoras huecas bajo la máquina.
  - Posicionadores de las fuentes y patas cuadradas con placas.
  - Gatos de elevación dobles con ojal.
  - Marcos de puerto y soportes de blanco, faldón de lamas, anillo de tuberías y chaflán de 45 mm en las aristas del yugo.
  - Autoblindaje opcional con las piezas de 4.645 × 2.770 × 2.575 m de la ficha.
- **Validación contra la ficha** (sin ajustar nada a estos números):
  - Envolvente de 1.92 × 1.92 × 1.83 m frente a 1.9 × 1.9 × 1.8 m.
  - Masa de 18.3 t frente a 18 t.
- **Blancos Nirta F-500:** cavidad de 2.9 ml según la ficha (antes 1.3 ml estimada).
- **Vacío y Molflow+:**
  - Superficies del vacío clasificadas por material vecino; las 4 bombas suman 1774 cm².
  - STL y parámetros para MolflowCLI, y lector de sus resultados.
  - Pérdidas de H⁻ por gas residual calibradas a la transmisión medida, más puntos de impacto de H⁰ en la pared.
  - *Stripping* de Lorentz verificado despreciable: 4·10⁻¹⁰ /m a 18 MeV.
- **Vacío MC:** ahora incluye las perforaciones de los colimadores y el interior de las bombas y de las carcasas RF.
- **Render Cycles 4K/8K** (`blender/render_kiube.py`): vista exterior tipo ficha y vista abierta, con la mitad superior elevada.

## Qué cambió en v0.4

- **Impresión 3D vía Blender:**
  - `kiube/print3d.py` agrupa las piezas en 41 cuerpos imprimibles y escala, valida y exporta STL, `.glb` y manifest.
  - `blender/blender_prep_impresion.py` agrega pines de registro, orienta, divide por cama, acomoda en placas y exporta.
  - Probado con Blender 5.0.1: STL estancos, verificados en Blender y con trimesh.
- **Búnker de los modelos MC:** pasa a las medidas del plano BIND de Radiopharma-OS (`plano_bunker.py`). La sala es de 5.9×4.7 m (estimada sobre plano) y los muros de 2.11 m (medido el oeste). El laberinto aún no está modelado.

## Qué cambió en v0.3

- **Modelos Monte Carlo** (`--lod mc`):
  - Salen DAGMC para OpenMC y DAG-MCNP y GDML para Geant4, desde una sola malla Gmsh conforme.
  - Validación en Geant4: 0 solapes y 0 tracks atascados.
  - La fuente de neutrones usa el espectro H₂¹⁸O de Bonvin Fig. 7, extraído de los vectores del PDF (4.93·10⁻³ n/p).
- **Stripper:** su borde interior queda en R_ex − ΔR/2, con el criterio del TR24 (Kazarinov 2019). ΔE = 98.5 keV/vuelta, ΔR = 1.30 mm y ~183 vueltas.
- **Corrección geométrica:** la ranura del inserto era tangente al tapón F. Eso creaba una arista cuchillo no-manifold, que rompe Geant4 y DAGMC; ahora penetra 5 mm.
- **Material nuevo:** concreto Portland (NIST) con trazas de Eu y Co.

## Qué cambió en v0.2

| Elemento | v0.1 | v0.2 | Fuente |
|---|---|---|---|
| Arquitectura | sector = cuña bajo cada polo, valle profundo 580 mm | **A** tapa · **B** yugo lateral · **C** disco en 4 piezas · **D** polo delgado · **E** inserto · **F** tapón | T04F1 |
| Gap de valle | 580 mm | **132 mm** (B_v/B_h = 0.21) | S09 + T04F1 |
| Campo colina / valle | 1.95 / 0.77 T (EST) | **2.15 / 0.45 T** | S09 (mapa medido) |
| RF | 40.65 MHz | 40.646 MHz; f_rev del modelo = 20.328 MHz (objetivo 20.323) | T04F8 |
| Ancho efectivo de colina | = geométrico | **46.7°** resuelto por isocronía | CALC |
| Tapa A | 140 mm | 166 mm | T04F1 |
| Radio interior yugo B | 745 mm | **778 mm** | T04F2 |
| Chaflán de esquinas | 330 mm | **263 mm** (las 4) | T04F2 |
| Polo | 1 escalón + arco | escalón 1 R = 529 · escalón 2 L = 512, arco descentrado R = 328 (Rh/Lh = 0.64) · lados desplazados 10 mm | T04F2/F3, PAT |
| Ranura del inserto | 70 mm | **52 mm**, del tapón al borde | T04F2 |
| Agujeros de bombeo | r 95 @ 300 | **r 127 @ 375**, atraviesan A y C | T04F2 |
| Bobinas | 590–730 × 160 | **568–759** × 190, apoyadas en A | T04F2 |
| Cámara | 552–577 | 548–566 | T04F2 |
| Puertos de extracción | 4 | **8** (2 por cara, a +17.2° y +44.5° del centro de colina) | T04F2 |
| Strippers | azimut supuesto | **ajustados por trazado** para que el haz pase por los puertos medidos: −14.2° y +4.9° | CALC |
| Puertos de bobina | 1 genérico | a 0° y 180° | T04F2 |
| Brazos de stripper | radiales (chocaban con el colimador vecino) | inclinados −25° | TUD03 Fig.3 |

### Validaciones cruzadas

- **Strippers:** el ajuste a los puertos medidos da −14.2° y +4.9°. En TUP03 Fig. 2 los strippers del KIUBE aparecen a unos −10° y +4°.
- **Cortes correctores:** el cruce del borde del 2.º sistema cae a **+17.8°**. TUP04 Fig. 2 muestra el corte a +17°. El 1.º sistema sale a −2.3° contra −14° medido, la principal discrepancia abierta.
- **Masas:** hierro del imán 13.5 t y bobinas 1.19 t c/u (18/9: 20 t de hierro, 1.22 t por bobina).

## Pendiente, por impacto

1. **Imán:** cerrado en v0.21 con 19 kW (IPAC26) y las cotas de las patentes. Queda la discrepancia con «< 11 kW» de TUD03T, registrada. **Curva B-H medida y altura de la bobina.** Calibrado en v0.18: 116 kA·v frente a 105 (potencia). Lo que falta es no único: J_s (2.12 T, CALC), la tapa (181, en el borde de la tolerancia) y la altura de la bobina (190, EST) se compensan entre sí. Una curva del lote o una foto con escala de la bobina lo cierran.
2. **Blindaje del búnker del proyecto.** Emisión hacia adelante de los blancos, búnker real (5.9 × 4.7 m, laberinto sur, techo) y humedad del hormigón. **Necesita el plano completo.** El benchmark contra el SRG ya está hecho (v0.13).
3. **RF.** Calibrada a 40.646 MHz, con puente central, sintonizador y acoplador (v0.19). Su geometría es EST y la calibración no es única (altura del dee, Ø del stem, gap del dummy dee, ancho del puente). Falta:
   - confirmar si IBA une los dees en el centro o fuera (KEY22, Fig. 5: barra entre stems);
   - confirmar el tipo de acoplador (lazo o capacitivo);
   - medir el factor de Q real (0.65, EST).

   Un plano o una foto del interior de la carcasa RF lo cerrarían.
4. **Cotas estimadas restantes (EST):** 36, listadas por impacto en `docs/COTAS_EST.md` (v0.26). Las de más impacto son la altura de la bobina, la abertura del dee, el factor de Q real, la velocidad de bombeo, el flujo de H₂, el piso del búnker y las cotas del blanco Nirta. Se cierran con el System Description de IBA (ref. 58027), la ficha del Nirta o fotos con escala.
5. **Valle (100 mm, Gv 200 de la patente), tapón (−15 mm), tapa (181 mm) y cuña del corrector (78 × 39 mm, patente):** ajustados con FE contra KEY22 y S09 (v0.21); falta confirmarlos con un plano.
6. **Región central:** puntas, cuello y posición de las fuentes tomados del mapa Opera de IBA (v0.22). Gap fuente–puller (2.8 mm, v0.23) y forma del puller (10 × 7 × 3 mm, v0.24) optimizados, y puntas a ±8 mm. Rendija de la chimenea elegida por corriente útil (6 mm, v0.25). Faltan la forma de las narices, y saber si el mapa está espejado (orientación azimutal de las chimeneas).
7. **Vacío:** desgasificación calibrada con el bombeo medido (v0.20). Faltan la velocidad real de las 4 VHS-6 para H₂ tras válvula y baffle, y el flujo de H₂ de la fuente. Falta repetir la comparación con Molflow (hay que volver a aportar `molflowCLI`).
8. **Motores comerciales:** Opera-3D/TOSCA (paquete en `out/opera3d/`, con la curva y el NI de v0.18). Si hay licencia, faltan los scripts `.comi` de la versión concreta.
