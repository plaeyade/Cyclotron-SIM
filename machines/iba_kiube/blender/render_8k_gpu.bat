@echo off
REM Render 8K del Cyclone KIUBE en GPU NVIDIA (OptiX). Doble clic desde el Explorador.
REM Busca Blender 4.4+ en Program Files; o define BLENDER antes de ejecutar.
cd /d "%~dp0.."
if not defined BLENDER (
  for /d %%D in ("C:\Program Files\Blender Foundation\Blender *") do if exist "%%D\blender.exe" set BLENDER="%%D\blender.exe"
)
if not defined BLENDER (
  where blender >nul 2>nul && set BLENDER=blender
)
if not defined BLENDER (
  echo No encontre blender.exe. Define BLENDER=ruta\blender.exe y vuelve a ejecutar.
  pause
  exit /b 1
)
echo Usando %BLENDER%
if not exist out\render\final mkdir out\render\final
%BLENDER% -b -P blender\render_kiube.py -- --escena out\render\scene.json --salida out\render\final\kiube_exterior_8k.png --res 7680x4320 --muestras 128 --gpu auto
%BLENDER% -b -P blender\render_kiube.py -- --escena out\render\scene.json --salida out\render\final\kiube_abierta_8k.png --res 7680x4320 --muestras 128 --gpu auto --vista abierta
if exist out\central_region\c18_9\escena\escena.json (
  %BLENDER% -b -P blender\render_region_central.py -- --escena out\central_region\c18_9\escena\escena.json --salida out\render\final\region_central_c18_9_4k.png --res 3840x2160 --muestras 256 --gpu auto
)
echo Listo: out\render\final\ (exterior 8K, abierta 8K, region central 4K)
pause
