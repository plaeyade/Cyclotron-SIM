#!/bin/bash
# Corre el benchmark SRG: espesores de muro × 2 semillas en paralelo (2 núcleos). Uso: correr_srg.sh N f "t1 t2 ..." [carpeta]
N=${1:-6000}; F=${2:-3.0}; TS=${3:-"100 150 200 220 230"}; D=${4:-out/mc_srg}
G=/home/claude/mamba/g4; export MAMBA_ROOT_PREFIX=/home/claude/mamba/root
cd /home/claude/kiube/$D && mkdir -p res
for t in $TS; do
  for s in 1 2; do
    /home/claude/mamba/bin/micromamba run -p $G /home/claude/kiube/mc_benchmark/srg_bench srg_config.txt kiube_srg.gdml \
      $t $N $((1000*t+s)) res/t${t}_s${s}.txt $F > res/t${t}_s${s}.log 2>&1 &
  done
  wait
  echo "hecho t=$t $(date +%H:%M:%S)"
done
echo FIN
