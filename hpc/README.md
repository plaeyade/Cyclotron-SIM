# Correr el Monte Carlo del búnker en una PC con muchos núcleos (p. ej. Alienware)

Cada lote de `bunker_iba` es un proceso de un solo hilo, independiente, con su semilla. El script lanza tantos lotes en paralelo como núcleos tengas y los reparte entre las caras; los resultados se suman después con `analizar_f90.py`. **La GPU no ayuda:** Geant4 transporta en CPU. Lo que acelera es el número de núcleos lógicos (y RAM: ~0.5–0.9 GB por proceso).

## Windows (Docker Desktop + WSL2)

```powershell
git clone https://github.com/plaeyade/Cyclotron-SIM.git
cd Cyclotron-SIM
docker build -t cyclotron-sim-mc -f hpc/Dockerfile .
mkdir salida
docker run -d --name mc_f90 --cpus 22 -v "${PWD}\salida:/salida" cyclotron-sim-mc todas 22 200 22   # 22 núcleos, 200 historias/lote, 22 lotes/cara
```

Variantes: `... cyclotron-sim-mc techo 0 200 32` (solo techo, 32 lotes); `... todas 0 200 0 td250bajoH.txt` (TD-250 bajo en hidrógeno); `... todas 0 200 0 td250std.txt`.

Se puede detener (Ctrl+C) y relanzar: los lotes terminados se saltan.

## Subir los resultados

Copia la carpeta `salida` (o solo los `*_regiones.txt` y `*_malla.bin`) y compártela; o súbela a una rama `resultados/<fecha>` de este repo. Con eso se corre `machines/iba_kiube/mc_benchmark/analizar_f90.py`.
