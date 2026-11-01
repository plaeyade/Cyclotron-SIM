"""F63 — Comparación cuantitativa CSG vs DAGMC (misma fuente, material y tally).

Replica el caso F62 en geometría CSG pura con la API Python de OpenMC:
- Misma cáscara: caja exterior 971×1549×400 cm, cavidad 211..760 / 116..1433 / 116..284 cm.
- Mismo hormigón Portland ordinario (composición F62, rho=2.3 g/cm3).
- Misma fuente puntual isotrópica en (485.5, 774.5, 200) cm con espectro BIND
  nspec_o18.tot (tabular, 163 bins).
- Mismo tally: flujo (track-length) en el volumen del muro.
Resultado de referencia DAGMC (F62): 595.88 ± 0.13 cm/partícula.
"""
from __future__ import annotations
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos

import json
import sys
from pathlib import Path

import numpy as np

E_DEF = Path(f"{_RPOS_RAIZ}/work/evidencia_v072/fase63")


def correr_csg(n_particles: int = 200_000, batches: int = 20,
               salida: Path = E_DEF) -> dict:
    import openmc

    spec = np.loadtxt(f"{_RPOS_RAIZ}/work/evidencia_v072/fase32/nspec_o18.tot",
                      comments="#")
    Ee = spec[:, 0]; w = spec[:, 1]
    m = w > 0
    Ee = Ee[m] * 1e6; w = w[m]

    mat = openmc.Material(name="Concreto")
    mat.set_density("g/cm3", 2.3)
    for el, wo in [("H", 0.010), ("O", 0.532), ("Na", 0.029), ("Al", 0.034),
                   ("Si", 0.337), ("Ca", 0.044), ("Fe", 0.014)]:
        mat.add_element(el, wo)

    # cáscara CSG idéntica a la DAGMC (cm)
    ox0, ox1 = openmc.XPlane(0.0), openmc.XPlane(971.0)
    oy0, oy1 = openmc.YPlane(0.0), openmc.YPlane(1549.0)
    oz0, oz1 = openmc.ZPlane(0.0), openmc.ZPlane(400.0)
    ix0, ix1 = openmc.XPlane(211.0), openmc.XPlane(760.0)
    iy0, iy1 = openmc.YPlane(116.0), openmc.YPlane(1433.0)
    iz0, iz1 = openmc.ZPlane(116.0), openmc.ZPlane(284.0)
    for s in (ox0, oy0, oz0):
        s.boundary_type = "vacuum"
    for s in (ox1, oy1, oz1):
        s.boundary_type = "vacuum"
    for s in (ix0, iy0, iz0, ix1, iy1, iz1):
        s.boundary_type = "transmission"
    outer = +ox0 & -ox1 & +oy0 & -oy1 & +oz0 & -oz1
    inner = +ix0 & -ix1 & +iy0 & -iy1 & +iz0 & -iz1
    muro = openmc.Cell(name="muro", region=outer & ~inner, fill=mat)
    cav = openmc.Cell(name="cavidad", region=inner)
    geom = openmc.Geometry([muro, cav])

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((485.5, 774.5, 200.0))
    src.angle = openmc.stats.Isotropic()
    src.energy = openmc.stats.Tabular(Ee, w, interpolation="linear-linear")

    st = openmc.Settings()
    st.run_mode = "fixed source"
    st.particles = n_particles
    st.batches = batches
    st.source = src

    t = openmc.Tally(name="flujo_muro")
    t.filters = [openmc.CellFilter([muro])]
    t.scores = ["flux"]
    tallies = openmc.Tallies([t])

    model = openmc.Model(geometry=geom, materials=[mat], settings=st,
                         tallies=tallies)
    import os, tempfile
    with tempfile.TemporaryDirectory() as td:
        os.environ["OPENMC_CROSS_SECTIONS"] = (
            f"{_RPOS_RAIZ}/work/sci_runtime/openmc_data/endfb80/"
            "cross_sections_sab.xml")
        sp = model.run(cwd=td)
        with openmc.StatePoint(sp) as s:
            tal = s.tallies[t.id]
            media = float(tal.mean.flatten()[0])
            sigma = float(tal.std_dev.flatten()[0])

    ref = 595.88; ref_s = 0.13
    res = {
        "fase": "F63_comparacion_csg_vs_dagmc",
        "estado": "completada_real",
        "csg": {"flujo_muro": media, "sigma": sigma,
                "particles": n_particles, "batches": batches},
        "dagmc_f62": {"flujo_muro": ref, "sigma": ref_s,
                      "particles": 200000, "batches": 20},
        "diferencia_relativa": abs(media - ref) / ref,
        "sigma_combinada": (sigma**2 + ref_s**2) ** 0.5,
        "mismo_modelo": "fuente puntual BIND + hormigon Portland + tally flux",
    }
    res["en_sigmas"] = abs(media - ref) / res["sigma_combinada"]
    Path(salida).mkdir(parents=True, exist_ok=True)
    (Path(salida) / "comparacion_csg_dagmc_f63.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    return res


if __name__ == "__main__":
    r = correr_csg()
    print(f"CSG: {r['csg']['flujo_muro']:.2f} ± {r['csg']['sigma']:.2f}")
    print(f"DAGMC: {r['dagmc_f62']['flujo_muro']:.2f} ± {r['dagmc_f62']['sigma']:.2f}")
    print(f"diferencia: {r['diferencia_relativa']*100:.2f} % "
          f"({r['en_sigmas']:.1f} σ)")
