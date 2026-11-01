# Cyclotron-SIM — imagen de cómputo Monte Carlo (Geant4 11.4.3 QGSP_BIC_HP + datos HP)
# Construir:  docker build -t cyclotron-sim-mc -f hpc/Dockerfile .
# Correr:     docker run -d --name mc_f90 --cpus 22 -v "${PWD}\salida:/salida" cyclotron-sim-mc todas 22 200 22
FROM mambaorg/micromamba:1.5.10
USER root
RUN micromamba install -y -n base -c conda-forge geant4=11.4.3 gxx cxx-compiler xerces-c python=3.12 numpy && micromamba clean -a -y
WORKDIR /opt/cyclotron-sim
COPY machines/iba_kiube/mc_benchmark /opt/cyclotron-sim/mc_benchmark
COPY machines/iba_kiube/resultados/mc_f90 /opt/cyclotron-sim/mc_f90
COPY hpc/correr_paralelo.py /opt/cyclotron-sim/correr_paralelo.py
RUN micromamba run -n base bash -c 'set -e; cd mc_benchmark; CXXB=$(command -v g++ || command -v x86_64-conda-linux-gnu-c++); echo "compilador: $CXXB"; \
    $CXXB -O3 -march=x86-64-v3 -std=c++17 bunker_iba.cc -o bunker_iba $(geant4-config --cflags) $(geant4-config --libs) -lxerces-c -Wl,-rpath,$CONDA_PREFIX/lib; \
    cd ../mc_f90; python -c "import lzma,glob,shutil; [shutil.copyfileobj(lzma.open(f), open(f[:-3],\"wb\")) for f in glob.glob(\"gdml/*.xz\")]"; ls -la gdml'
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/cyclotron-sim/correr_paralelo.py"]
