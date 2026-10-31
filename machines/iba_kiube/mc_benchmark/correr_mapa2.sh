#!/bin/bash
# v2 (importancia por profundidad en muros, techo y losa + malla de flujo 3 grupos): bloques de 2000 historias
# correr_lote.sh semilla0 paso presupuesto_s [f]
export MAMBA_ROOT_PREFIX=/home/claude/mamba/root
cd /home/claude/kiube/out/mc_f70v035/h2500_t1450_v2
s=$1; T0=$(date +%s); F=${4:-2.15}
while [ $(( $(date +%s)-T0 )) -lt $3 ]; do
  t0=$(date +%s)
  timeout $(( $3 - ($(date +%s)-T0) + 60 )) /home/claude/mamba/bin/micromamba run -p /home/claude/mamba/g4 /home/claude/kiube/mc_benchmark/bunker_mapa2 mapa_config.txt kiube.gdml 2000 $s res/r$s $F > res/r$s.log 2>&1
  echo "r$s rc=$? $(( $(date +%s)-t0 )) s"
  s=$(( s+$2 ))
done
