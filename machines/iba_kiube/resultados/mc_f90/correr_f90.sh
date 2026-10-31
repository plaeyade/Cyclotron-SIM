#!/bin/bash
# uso: correr_f90.sh config gdml cara N semilla1 [semilla2 ...]
export MAMBA_ROOT_PREFIX=/home/claude/mamba/root; eval "$(/home/claude/mamba/bin/micromamba shell hook -s bash)"; micromamba activate /home/claude/mamba/g4
cd /home/claude/kiube/out/mc_f90; B=/home/claude/kiube/mc_benchmark/bunker_iba
cfg=$1; g=$2; cara=$3; N=$4; shift 4
case $cara in este) a="3.0 1 1 2000 25";; oeste) a="3.0 1 -1 2000 27";; sur) a="3.0 0 1 1970 28";; techo) a="3.0 2 1 1380 26";; esac
tag=$(basename $cfg .txt)
for s in "$@"; do pre=res/${tag}_${cara}_s$s
  [ -s ${pre}_regiones.txt ] && continue
  ( $B $cfg $g $N $s ${pre}_tmp $a > ${pre}.log 2>&1 && mv ${pre}_tmp_malla.bin ${pre}_malla.bin && mv ${pre}_tmp_regiones.txt ${pre}_regiones.txt && rm -f ${pre}.log ) &
done; wait
