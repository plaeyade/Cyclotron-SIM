#!/bin/bash
# Corre un caso (espesor, semilla) del benchmark SRG. Uso: correr_caso.sh t semilla [N] [f] [carpeta]
t=$1; s=$2; N=${3:-6000}; F=${4:-3.0}; D=${5:-out/mc_srg}
G=/home/claude/mamba/g4; export MAMBA_ROOT_PREFIX=/home/claude/mamba/root
cd /home/claude/kiube/$D && mkdir -p res
/home/claude/mamba/bin/micromamba run -p $G /home/claude/kiube/mc_benchmark/srg_bench srg_config.txt kiube_srg.gdml \
  $t $N $((1000*t+s)) res/t${t}_s${s}.txt $F > res/t${t}_s${s}.log 2>&1
echo "hecho t=$t s=$s rc=$? $(date +%H:%M:%S)"
