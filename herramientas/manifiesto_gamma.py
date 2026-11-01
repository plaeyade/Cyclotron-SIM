"""Genera gamma-ingesta.json (contrato gamma/ingesta@1, ADR-0034 de GAMMA-AI) para machines/*/fuentes/.

Uso: python herramientas/manifiesto_gamma.py   (desde la raíz del repo)
Verificar con GAMMA:  gamma ecosistema validar gamma-ingesta.json --raiz-productor .
Sólo stdlib. Los metadatos por carpeta están en CATALOGO; un archivo sin entrada falla.
"""
import datetime
import hashlib
import json
import mimetypes
import pathlib
import subprocess

RAIZ = pathlib.Path(__file__).resolve().parents[1]

# Documentos de fabricante (IBA) → restricted; literatura técnica publicada → public.
FABRICANTE = ("IBA_", "MID_", "cbr_cyclone", "2024-10-16-Feedback")


def metadatos(p):
    maquina = p.parts[-3]
    fab = p.name.startswith(FABRICANTE)
    return {
        "dominio": f"ciclotron/{maquina}",
        "gobernanza": "restricted" if fab else "public",
        "emisor": "IBA" if fab else "literatura tecnica (congresos/revistas)",
        "licencia": "documentacion de fabricante; uso interno" if fab else "publicacion de acceso abierto; cita obligatoria",
        "idioma": "es" if "Propuesta" in p.name else "en",
    }


def main():
    docs = []
    for p in sorted(RAIZ.glob("machines/*/fuentes/*")):
        if not p.is_file():
            continue
        md = metadatos(p)
        data = p.read_bytes()
        rel = p.relative_to(RAIZ).as_posix()
        docs.append({
            "id_local": rel.rsplit(".", 1)[0],
            "ruta": rel,
            "sha256": hashlib.sha256(data).hexdigest(),
            "bytes": len(data),
            "mime": mimetypes.guess_type(p.name)[0] or "application/octet-stream",
            "titulo": p.stem.replace("_", " "),
            "dominio": md["dominio"],
            "gobernanza": md["gobernanza"],
            "tipo": "document",
            "emisor": md["emisor"],
            "licencia": md["licencia"],
            "idioma": md["idioma"],
        })
    try:
        commit = subprocess.run(["git", "rev-parse", "--short=12", "HEAD"], cwd=RAIZ,
                                capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        commit = None
    man = {"contrato": "gamma/ingesta@1",
           "productor": {"repo": "cyclotron-sim", **({"commit": commit} if commit else {})},
           "generado_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
           "documentos": docs}
    (RAIZ / "gamma-ingesta.json").write_text(json.dumps(man, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"gamma-ingesta.json: {len(docs)} documentos")


if __name__ == "__main__":
    main()
