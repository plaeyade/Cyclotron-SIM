#!/bin/bash
# Render 8K del Cyclone KIUBE en GPU (OptiX/CUDA). Uso: BLENDER=/ruta/blender bash blender/render_8k_gpu.sh
cd "$(dirname "$0")/.."
B=${BLENDER:-blender}
for V in exterior abierta; do
  "$B" -b -P blender/render_kiube.py -- --escena out/render/scene.json \
      --salida out/render/final/kiube_${V}_8k.png --res 7680x4320 --muestras 128 --gpu auto --vista $V
done
[ -f out/central_region/c18_9/escena/escena.json ] && "$B" -b -P blender/render_region_central.py -- \
    --escena out/central_region/c18_9/escena/escena.json --salida out/render/final/region_central_c18_9_4k.png \
    --res 3840x2160 --muestras 256 --gpu auto
