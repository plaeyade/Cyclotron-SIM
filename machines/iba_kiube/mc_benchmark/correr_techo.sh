#!/bin/bash
# Corre un caso del techo. Uso: correr_techo.sh carpeta t_techo_cm semilla [N] [f] [rho_techo]
D=$1; t=$2; s=$3; N=${4:-6000}; F=${5:-3.0}; R=${6:-}
G=/home/claude/mamba/g4; export MAMBA_ROOT_PREFIX=/home/claude/mamba/root
cd /home/claude/kiube/$D && mkdir -p res
tag=t${t}_s${s}${R:+_rho$R}
/home/claude/mamba/bin/micromamba run -p $G /home/claude/kiube/mc_benchmark/bunker_techo techo_config.txt kiube.gdml \
  $t $N $((1000*t+s)) res/$tag.txt $F $R > res/$tag.log 2>&1
echo "hecho $D $tag rc=$? $(date +%H:%M:%S)"
