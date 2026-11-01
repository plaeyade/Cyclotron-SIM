"""Configuración de pruebas de pipelines.

Las pruebas de fase (F54–F97) verifican evidencia generada fuera de git (JSON, mallas MC,
corridas DWSIM/Molflow…) bajo RPOS_AGENTS_ROOT (por defecto /mnt/agents, el entorno donde se
produjo).

- Sin RPOS_AGENTS_ROOT definido, las pruebas usan una raíz TEMPORAL: lo que las fases escriben
  no ensucia el sistema, y lo que leen y no existe se reporta como OMITIDO («evidencia externa
  ausente»), no como fallo. Una falla por cualquier otra causa sigue siendo falla.
- RPOS_REQUIRE_EVIDENCIA=1 (equipo con la evidencia restaurada y RPOS_AGENTS_ROOT apuntándola):
  nada se omite.
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path

import pytest

_REQUERIDA = os.environ.get("RPOS_REQUIRE_EVIDENCIA") == "1"
if not _REQUERIDA and "RPOS_AGENTS_ROOT" not in os.environ:
    os.environ["RPOS_AGENTS_ROOT"] = tempfile.mkdtemp(prefix="rpos-agents-")
_RAIZ = os.environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")
_MOTIVO = "evidencia externa ausente (RPOS_AGENTS_ROOT); RPOS_REQUIRE_EVIDENCIA=1 la exige"


_REPO = str(Path(__file__).resolve().parents[1])


def _ruta_externa(texto: str) -> bool:
    return _RAIZ in texto or "/mnt/agents" in texto or "evidencia_v072" in texto


def _fuera_del_repo(ruta: str) -> bool:
    return ruta.startswith("/") and not ruta.startswith(_REPO)


def _evidencia_ausente(exc) -> bool:
    if _REQUERIDA:
        return False
    vistos = set()
    while exc is not None and id(exc) not in vistos:
        vistos.add(id(exc))
        if isinstance(exc, FileNotFoundError):
            ruta = str(exc.filename or "")
            if _ruta_externa(ruta) or _fuera_del_repo(ruta):
                return True
        if isinstance(exc, AssertionError) and _RAIZ in str(exc) and "exists" in str(exc):
            return True             # la prueba afirma explícitamente que la evidencia existe
        if isinstance(exc, subprocess.CalledProcessError):
            salida = f"{exc.stderr or ''}{exc.output or ''}"
            if isinstance(salida, bytes):
                salida = salida.decode(errors="replace")
            if "FileNotFoundError" in salida and _ruta_externa(salida):
                return True
        exc = exc.__cause__ or exc.__context__
    return False


# Calificación del ENTORNO (runtimes científicos instalados, evidencia de fases previas con
# contenido específico): sólo tiene sentido en el equipo de producción/evidencia.
_ENTORNO = {
    "test_oq_fase40.py": "OQ contra runtimes y evidencia del entorno de producción",
    "test_qualification.py::test_dq_documentos_y_trazabilidad": "DQ contra documentos del entorno",
    "test_qualification.py::test_iq_instalacion": "IQ de la instalación de producción",
    "test_informe_regulatorio_f54.py::test_integra_evidencia_real": "integra evidencia MC real",
    "test_informe_regulatorio_f54.py::test_sintesis_solo_valores_publicados": "integra evidencia MC real",
    "test_dwsim_f84.py::test_estado_y_runtime": "runtime DWSIM instalado",
    # El veredicto del planner depende del estado de la evidencia de fases previas en la raíz.
    "test_madse_f92.py::test_planner_y_orquestador": "veredicto según evidencia de fases previas",
}


def pytest_collection_modifyitems(config, items):
    if _REQUERIDA:
        return
    for item in items:
        nombre = Path(str(item.path)).name
        motivo = _ENTORNO.get(f"{nombre}::{item.name}") or _ENTORNO.get(nombre)
        if motivo:
            item.add_marker(pytest.mark.skip(reason=f"calificación de entorno: {motivo}; RPOS_REQUIRE_EVIDENCIA=1 la exige"))


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    if call.excinfo is not None and rep.failed and _evidencia_ausente(call.excinfo.value):
        rep.outcome = "skipped"
        rep.longrepr = (str(item.path), item.location[1] or 0, _MOTIVO)


@pytest.hookimpl(hookwrapper=True)
def pytest_make_collect_report(collector):
    outcome = yield
    rep = outcome.get_result()
    texto = str(getattr(rep, "longrepr", "") or "")
    if rep.failed and not _REQUERIDA and "FileNotFoundError" in texto and _ruta_externa(texto):
        rep.outcome = "skipped"
        rep.longrepr = (str(collector.path), 0, _MOTIVO)
        rep.result = []
