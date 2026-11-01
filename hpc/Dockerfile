# Cyclotron-SIM — imagen de cómputo Monte Carlo (Geant4 11.4.3 QGSP_BIC_HP + datos HP)
# Construir:  docker build -t cyclotron-sim-mc -f hpc/Dockerfile .
# Correr:     docker run --rm -v "%cd%/salida:/salida" cyclotron-sim-mc todas 0      (todas las caras, todos los núcleos)
FROM mambaorg/micromamba:1.5.10
USER root
RUN micromamba install -y -n base -c conda-forge geant4=11.4.3 cxx-compiler xerces-c python=3.12 numpy xz && micromamba clean -a -y
SHELL ["/usr/local/bin/_dockerfile_shell.sh"]
WORKDIR /opt/cyclotron-sim
COPY machines/iba_kiube/mc_benchmark /opt/cyclotron-sim/mc_benchmark
COPY machines/iba_kiube/resultados/mc_f90 /opt/cyclotron-sim/mc_f90
COPY hpc/correr_paralelo.py /opt/cyclotron-sim/correr_paralelo.py
RUN cd mc_benchmark && g++ -O3 -march=x86-64-v3 -std=c++17 bunker_iba.cc -o bunker_iba $(geant4-config --cflags) $(geant4-config --libs) -lxerces-c \
 && cd ../mc_f90 && for f in gdml/*.xz; do xz -dk "$f"; done
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/cyclotron-sim/correr_paralelo.py"]
