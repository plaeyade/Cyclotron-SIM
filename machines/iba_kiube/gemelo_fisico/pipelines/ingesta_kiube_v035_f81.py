#!/usr/bin/env python3
"""F81 — Ingesta multiformato del gemelo kiube v0.35 + verificación cruzada contra KIUBE-V04, F75, SRG y KIUBE-OM.

Entrada: paquete generado en el repositorio kiube (herramientas/paquete_rpos.py → out/rpos_v035/):
  manifiesto.json (geometría STEP/h5m/GDML/STL/3MF/GLB/Molflow, mapas npz, curvas JSON, parámetros, telemetría simulada;
  cada artefacto con md5 y metadatos), verificacion_cruzada.json, registro_reetiquetado.jsonl, impacto_cadena.json,
  telemetria_simulada/.
Esta fase:
  1. verifica el md5 de CADA artefacto contra el archivo (integridad; nada se acepta sin hash coincidente);
  2. incorpora la verificación cruzada SIN reescribir historia: los supuestos KIUBE-V04 no se sobrescriben; el
     re-etiquetado es un registro append-only (valor y etiqueta anteriores conservados). Solo sube a
     `fuente_modelo_oficial_IBA (vía usuario)` lo que coincide dentro de la tolerancia declarada Y tiene fuente IBA en
     el modelo; lo demás se reporta como discrepancia cuantificada (Δ, Δ/ref);
  3. marca como afectadas (append-only) las fases que dependen de la geometría de cámara (F66/F74/F67, F65) y del voxel
     (F70 y derivadas), con prioridad; no modifica sus evidencias;
  4. registra la telemetría simulada como `simulado_declarado`. No entra como medida al Calibration Layer: la puerta
     17.21 rechaza origen != "real" (memoria 3).
Reglas canónicas: simulacion_teorica, pq no_declarada_nunca (14.9), cifras con fuente (49.11), append-only ALCOA+,
recomendación != decisión operacional (18.10).
Rutas: KIUBE_RAIZ (repo kiube), RPOS_EVIDENCIA (evidencia_v072).
"""
_RPOS_RAIZ = __import__("os").environ.get("RPOS_AGENTS_ROOT", "/mnt/agents")  # raíz de evidencia/datos externos
import datetime, hashlib, json, os, pathlib

KIUBE = pathlib.Path(os.environ.get("KIUBE_RAIZ", f"{_RPOS_RAIZ}/work/kiube_v035/kiube"))
PAQ = KIUBE / "out/rpos_v035"
EV = pathlib.Path(os.environ.get("RPOS_EVIDENCIA", f"{_RPOS_RAIZ}/work/evidencia_v072"))
OUT = EV / "fase81"
ETQ_IBA = "fuente_modelo_oficial_IBA (vía usuario)"
FASES_CAMARA = {"F66", "F74", "F67", "F65"}
FASES_VOXEL = {"F70"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def append_unico(path, filas, clave):
    """Añade filas a un .jsonl sin tocar las existentes (dedup por clave)."""
    previas = set()
    if path.exists():
        previas = {json.loads(l)[clave] for l in path.read_text().splitlines() if l.strip()}
    n = 0
    with open(path, "a") as f:
        for r in filas:
            if r[clave] in previas:
                continue
            f.write(json.dumps(r, ensure_ascii=False, default=float) + "\n"); n += 1
    return n


def ejecutar():
    OUT.mkdir(parents=True, exist_ok=True)
    man = json.load(open(PAQ / "manifiesto.json"))
    # 1) integridad: md5 de cada artefacto
    malos, faltan, formatos = [], [], {}
    for a in man["artefactos"]:
        p = KIUBE / a["ruta"]
        formatos[a["formato"]] = formatos.get(a["formato"], 0) + 1
        if not p.exists():
            faltan.append(a["ruta"]); continue
        if md5(p) != a["md5"]:
            malos.append(a["ruta"])
    # 2) verificación cruzada y re-etiquetado
    ver = json.load(open(PAQ / "verificacion_cruzada.json"))
    filas = ver["filas"]
    reg_k = [json.loads(l) for l in (PAQ / "registro_reetiquetado.jsonl").read_text().splitlines() if l.strip()]
    for r in reg_k:                                  # coherencia: ningún re-etiquetado sin coincidencia + fuente IBA
        f = next(x for x in filas if x["id"] == r["id"])
        if r["accion"] == "re_etiquetado":
            assert f["veredicto"] == "coincide" and f["fuente_IBA"], r["id"]
            assert r["etiqueta_nueva"] == ETQ_IBA and r["etiqueta_anterior"].startswith("supuesto_declarado")
    n_reg = append_unico(OUT / "registro_reetiquetado.jsonl",
                         [dict(r, fase=81, origen_paquete=man["version"]) for r in reg_k], "clave")
    # 3) impacto en cadena → marcas append-only
    imp = json.load(open(PAQ / "impacto_cadena.json"))["afectadas"]
    hoy = datetime.datetime.now().isoformat(timespec="seconds")
    marcas = []
    for x in imp:
        k = hashlib.md5(json.dumps([x["fase"], x["motivo"]], ensure_ascii=False).encode()).hexdigest()
        marcas.append(dict(clave=k, fecha=hoy, fase_afectada=x["fase"], prioridad=x["prioridad"], motivo=x["motivo"],
                           accion_recomendada=x["accion"], estado="afectada_pendiente_re_corrida",
                           nota="recomendación (18.10); la evidencia original de la fase no se modifica"))
    n_mar = append_unico(OUT / "marcas_afectacion.jsonl", marcas, "clave")
    # 4) telemetría simulada declarada
    tel = json.load(open(PAQ / "telemetria_simulada/telemetria_simulada_meta.json"))
    csvs = sorted((PAQ / "telemetria_simulada").glob("*.csv"))
    origen_ok = all(l.rstrip().endswith("simulado_declarado") for c in csvs for l in c.read_text().splitlines()[1:] if l.strip())

    disc = [dict(id=f["id"], parametro=f["parametro"], referencia=f["valor_referencia"], modelo=f["valor_modelo"],
                 unidad=f["unidad"], delta=f["delta"], delta_rel=f["delta_rel"], nota=f["nota"])
            for f in filas if f["veredicto"] == "discrepancia"]
    out = {
        "fase": 81, "version": "3.7.0", "tipo": "ingesta_verificacion_cruzada", "pq": "no_declarada_nunca",
        "fecha": hoy, "titulo": "Ingesta multiformato kiube v0.35 + verificación cruzada V04/F75/SRG/KIUBE-OM",
        "paquete": {"modelo": man["modelo"], "version": man["version"], "md5_manifiesto": md5(PAQ / "manifiesto.json"),
                    "n_artefactos": man["n_artefactos"], "por_categoria": man["por_categoria"], "por_formato": formatos,
                    "bytes_total": man["bytes_total"]},
        "integridad": {"verificados": len(man["artefactos"]) - len(malos) - len(faltan), "md5_distinto": malos,
                       "faltantes": faltan, "ok": not malos and not faltan},
        "verificacion": {"n": len(filas),
                         "coinciden": [f["id"] for f in filas if f["veredicto"] == "coincide"],
                         "discrepancias": disc,
                         "no_comparables": [dict(id=f["id"], nota=f["nota"]) for f in filas if f["veredicto"] == "no_comparable"],
                         "re_etiquetadas": [f["id"] for f in filas if f["re_etiquetar"]],
                         "regla": ver["regla_reetiquetado"]},
        "registro": {"nuevas_lineas_reetiquetado": n_reg, "nuevas_marcas_afectacion": n_mar},
        "impacto": {"afectadas": imp,
                    "camara_difiere": any(x["fase"] in FASES_CAMARA for x in imp),
                    "voxel_afectado": any(x["fase"] in FASES_VOXEL for x in imp)},
        "telemetria_simulada": {"origen": tel["origen"], "csv": [c.name for c in csvs], "origen_por_fila_ok": origen_ok,
                                "canales": tel["canales"], "calibration_layer": tel["calibration_layer"],
                                "discrepancia_conocida": tel["discrepancia_conocida"]},
        "advertencias": [
            "La telemetría es del gemelo (simulado_declarado), no haz real: memoria 3 intacta; la puerta 17.21 del Calibration Layer rechaza origen != real.",
            "Re-etiquetar exige coincidencia Y fuente IBA en el modelo; una coincidencia con una cota EST/CALC no confirma nada.",
            "El canónico RPOS cita KIUBE-SRG rev. F (2019); el modelo kiube usa rev. G (2023). Las cifras comparadas coinciden en ambas revisiones (Tabla 1).",
            "KIUBE-OM se compara a través de su transcripción en docs/fuentes_oficiales_kiube.md (el documento no está en el repositorio kiube).",
        ],
    }
    json.dump(out, open(OUT / "ingesta_kiube_v035_f81.json", "w"), indent=1, ensure_ascii=False, default=float)
    with open(EV / "registro_fases.jsonl", "a") as f:
        f.write(json.dumps({"fase": 81, "version": "3.7.0", "fecha": hoy, "ok_integridad": out["integridad"]["ok"],
                            "discrepancias": len(disc), "afectadas": [x["fase"] for x in imp]}, ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    o = ejecutar()
    print(json.dumps({k: o[k] for k in ("paquete", "integridad", "registro")}, ensure_ascii=False, default=float)[:1500])
    print("discrepancias:", [(d["id"], d["delta_rel"]) for d in o["verificacion"]["discrepancias"]])
    print("afectadas:", [(x["fase"], x["prioridad"]) for x in o["impacto"]["afectadas"]])
