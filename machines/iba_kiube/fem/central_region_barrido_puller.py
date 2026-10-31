"""Barrido de la forma del puller de la región central del KIUBE (preset 'kiube', v0.24).

Puller = dos mordazas (arriba y abajo, unidas al dee) frente a la rendija de la fuente:
  puller_slit_h  abertura vertical entre mordazas [mm]
  puller_w       ancho tangencial [mm]
  puller_len     profundidad (cambia el alcance chimenea+gap+puller; se compensa src_dtheta para mantener la
                 fuente a 40.3° del borde de entrada, como en el mapa de IBA)
Cada caso: malla fina, 6 vueltas, N partículas (resultado completo en out/central_region/puller/<tag>/).
Uso: python fem/central_region_barrido_puller.py tag clave=valor [clave=valor ...]
"""
import json, math, pathlib, subprocess, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
ANG_IBA = math.degrees((3.0 + 2.8 + 3.0) / 12.5)          # 40.3° con el preset v0.23


def correr(tag, sets, n=400):
    sets = dict(sets)
    L = float(sets.get("puller_len", 3.0)); gap = float(sets.get("gap_sp", 2.8)); r = 12.5
    rc = float(sets.get("chimney_r", 3.0))                     # v0.32: el centro de la chimenea queda a 40.3° del borde
    dth = math.degrees((rc + gap + L) / r) - ANG_IBA
    fu = json.dumps([dict(especie="H-", h=2, dee=1, src_r=r, src_dtheta=round(dth, 3))])
    out = ROOT / "out/central_region/puller" / tag
    cmd = [sys.executable, str(ROOT / "fem/central_region.py"), "--preset", "kiube", "--n_haz", str(n), "--out", str(out),
           "--set", *[f"{k}={v}" for k, v in sets.items()], f"fuentes={fu}"]
    r_ = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if r_.returncode != 0:
        return dict(tag=tag, error=r_.stderr[-600:])
    R = json.load(open(out / "resultado.json")); e = R["especies"][0]
    nn = sum(e["perdidas_haz"].values()); p = e["perdidas_haz"]["vive"] / nn
    ok = [x for x in e["barrido"] if x["sobrevive"] and x["centro_mm"]]
    dm = sum(math.hypot(x["centro_mm"][0], x["centro_mm"][1]) for x in ok) / max(len(ok), 1)
    return dict(tag=tag, sets=sets, src_dtheta=round(dth, 3), transmision_pct=100 * p,
                error_pct=100 * math.sqrt(p * (1 - p) / nn), aceptancia_deg=e["ancho_aceptancia_deg"],
                descentrado_medio_mm=dm, perdidas={k: 100 * v / nn for k, v in e["perdidas_haz"].items()},
                guia_CAS=e.get("fraccion_ventana_en_guia_CAS"))


if __name__ == "__main__":
    tag = sys.argv[1]
    sets = dict(kv.split("=") for kv in sys.argv[2:])
    f = correr(tag, sets)
    db = ROOT / "out/central_region/puller/barrido.json"
    filas = json.load(open(db)) if db.exists() else []
    filas = [x for x in filas if x["tag"] != tag] + [f]
    json.dump(filas, open(db, "w"), indent=1, default=float)
    print("[puller]", json.dumps(f, default=float))
