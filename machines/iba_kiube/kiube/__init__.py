import yaml, pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
def load_params(path=None):
    raw = yaml.safe_load(open(path or ROOT / "params.yaml"))
    return {g: {k: v["v"] for k, v in d.items()} for g, d in raw.items()}, raw
