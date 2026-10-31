import json, pathlib, sys
import pytest
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "fem"))


@pytest.fixture(scope="session")
def PE():
    from kiube import load_params
    from kiube.beam import extraction
    P, raw = load_params()
    return P, extraction(P), raw


@pytest.fixture(scope="session")
def mats():
    import yaml
    return yaml.safe_load(open(ROOT / "materials.yaml"))


@pytest.fixture(scope="session")
def lods(PE):
    from kiube.geometry import build
    P, E, _ = PE
    return {l: build(P, E, l)[0] for l in ("print", "neutronics", "fem_struct")}
