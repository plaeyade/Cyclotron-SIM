#!/bin/bash
# Compila los benchmarks Geant4 (requiere Geant4 ≥ 11.2 con datos HP; probado con 11.4.3 vía micromamba).
set -e; cd "$(dirname "$0")"
for f in bunker_iba bunker_mapa bunker_mapa_v1b bunker_techo losa_tvl srg_bench; do
  [ -f $f.cc ] || continue
  g++ -O2 -std=c++17 $f.cc -o $f $(geant4-config --cflags) $(geant4-config --libs) -lxerces-c -Wl,-rpath,"$(geant4-config --prefix)/lib"
  echo "ok $f"
done
