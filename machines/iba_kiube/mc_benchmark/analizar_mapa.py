"""Análisis de las corridas de bunker_mapa (equivalente v0.35 de RPOS F70) → mapa H*(10), métricas y activación.

Combina corridas independientes (semillas) → media y desviación por vóxel (como F77). Normaliza a KIUBE 180 de
doble haz: 2 × 90 µA en T1 y T5 (PROP), 18O(p,n) con el rendimiento de BON21 (4.93e-3 n/p).
Activación de saturación (reacción de captura, convención de Høgdahl): R = N·(σ0·T_th + I0·T_epi/Δu), con T_th la
longitud de traza ponderada 1/v (E < 0.5 eV) y T_epi la de 0.5 eV–0.1 MeV (Δu = ln(1e5/0.5)).
Uso: python mc_benchmark/analizar_mapa.py out/mc_f70v035/h2500_t1450 [patrón=r*]
"""
from __future__ import annotations
import glob, json, math, pathlib, struct, sys
import numpy as np
import yaml
ROOT = pathlib.Path(__file__).resolve().parents[1]
QE = 1.602176634e-19; NA = 6.02214076e23

# secciones de captura (b): σ0 térmica 2200 m/s e integral de resonancia I0. Mughabghab, Atlas of Neutron Resonances
# (2006) para Co59, Cu63; Eu151→152g e I0 de Ag109→110m: orden de tabla, EST (pendiente_fuente de cita exacta)
XS = {"Co60": dict(el="Co", A=58.933, ab=1.0, s0=37.18, I0=75.9, t12_d=1925.28, fuente="Mughabghab 2006"),
      "Eu152": dict(el="Eu", A=151.96, ab=0.4781, s0=5900.0, I0=3300.0, t12_d=4943.0, fuente="EST (Eu151→152g, orden de tabla)"),
      "Cu64": dict(el="Cu", A=63.546, ab=0.6915, s0=4.50, I0=4.97, t12_d=0.5291, fuente="Mughabghab 2006"),
      "Ag110m": dict(el="Ag", A=107.87, ab=0.4816, s0=4.7, I0=70.0, t12_d=249.8, fuente="σ0 Mughabghab 2006; I0 EST")}
# qué región aporta cada isótopo (F72: Co60 en acero, Eu152 en hormigón, Ag110m y Cu64 en cobre)
CANALES = {"Co60": ["hormigon", "steel_1010", "stainless_304L"], "Eu152": ["hormigon"],
           "Cu64": ["copper_coil_homog", "copper_ofe", "copper_dee_caja"], "Ag110m": ["copper_coil_homog", "copper_ofe", "copper_dee_caja"]}
DU = math.log(1e5 / 0.5)


def leer_malla(f):
    b = open(f, "rb").read()
    nx, ny, nz = struct.unpack_from("<3i", b, 0)
    x0, y0, z0, h, N, _, _ = struct.unpack_from("<7d", b, 12)
    o = 12 + 56; nv = nx * ny * nz
    mn = np.frombuffer(b, "<f4", nv, o).reshape(nx, ny, nz)
    mg = np.frombuffer(b, "<f4", nv, o + 4 * nv).reshape(nx, ny, nz)
    d = dict(n=mn.astype(np.float64), g=mg.astype(np.float64), N=N, x0=x0, y0=y0, z0=z0, h=h, dims=(nx, ny, nz))
    if len(b) >= o + 20 * nv:                                       # v2: fluencia de neutrones en 3 grupos
        d["flux3g"] = np.stack([np.frombuffer(b, "<f4", nv, o + (2 + k) * 4 * nv).reshape(nx, ny, nz) for k in range(3)])
    return d


def leer_regiones(f):
    R = {}
    for l in open(f):
        if l.startswith("#"):
            continue
        t = l.split()
        if t[0] == "maquina":
            R[t[1]] = np.array(list(map(float, t[2:6])))
        else:
            k = ("hormigon", t[1], int(t[2]))
            R[k] = np.array(list(map(float, t[3:7])))
    return R


def main(d, patron="r*"):
    d = pathlib.Path(d).resolve()
    geo = json.load(open(d / "geometria.json"))
    cfg = {l.split()[0]: l.split()[1:] for l in open(d / "mapa_config.txt") if l.strip() and not l[0].isdigit() and not l.startswith("spectrum")}
    mallas = sorted(glob.glob(str(d / "res" / f"{patron}_malla.bin")))
    assert mallas, "sin corridas"
    M = [leer_malla(f) for f in mallas]
    V = M[0]["h"] ** 3 / 1000.0                                     # cm³ por vóxel
    # H*(10) por n fuente en cada vóxel [pSv]: suma(L·h)/N/V
    H = np.array([(m["n"] + m["g"]) / m["N"] / V for m in M]); Hn = np.array([m["n"] / m["N"] / V for m in M])
    Hm, Hs = H.mean(0), H.std(0, ddof=1) / math.sqrt(len(M)) if len(M) > 1 else np.zeros_like(H[0])
    Ntot = sum(m["N"] for m in M)
    # tasa de fuente: 2 × 90 µA (PROP) con el rendimiento BON21
    I_uA = 90.0; Y = geo["yield_BON21"]
    S = 2 * I_uA * 1e-6 / QE * Y                                     # n/s total (T1 + T5)
    k_uSvh = S * 3600 * 1e-6                                         # pSv por n → µSv/h
    nx, ny, nz = M[0]["dims"]; x0, y0, z0, h = M[0]["x0"], M[0]["y0"], M[0]["z0"], M[0]["h"]
    xc = x0 + (np.arange(nx) + 0.5) * h; yc = y0 + (np.arange(ny) + 0.5) * h; zc = z0 + (np.arange(nz) + 0.5) * h
    iz0 = int(np.argmin(np.abs(zc - 0.0)))                           # plano medio (altura de los blancos)
    ix_, iy_, ih_ = map(float, cfg["vault_inner"]); wall = float(cfg["wall"][0]); fz = float(cfg["floor_z"][0])
    slab = float(cfg["slab"][0]); roof = float(cfg["roof"][0])
    hx, hy = ix_ / 2, iy_ / 2
    X, Yg, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    dentro = (np.abs(X) < hx) & (np.abs(Yg) < hy) & (Z > fz) & (Z < fz + ih_)
    bloque = (np.abs(X) < hx + wall) & (np.abs(Yg) < hy + wall) & (Z > fz - slab) & (Z < fz + ih_ + roof)
    fuera = ~bloque

    def caja(c, semi_mm, campo=H):
        """media de H en un cubo centrado en c (semilado semi_mm) por corrida → media y error estándar entre corridas"""
        sel = (np.abs(xc - c[0]) <= semi_mm, np.abs(yc - c[1]) <= semi_mm, np.abs(zc - c[2]) <= semi_mm)
        v = np.array([Hk[np.ix_(*sel)].mean() for Hk in campo])
        return float(v.mean() * k_uSvh), (float(v.std(ddof=1) / math.sqrt(len(v)) * k_uSvh) if len(v) > 1 else None)

    T1 = geo["sources"][0]; p1 = np.array(T1["pos_mm"]); b = math.radians(T1["beam_dir_deg"])
    srg = {}
    for ang in (0, 90, 180):
        a = b + math.radians(ang); q = p1 + 1000 * np.array([math.cos(a), math.sin(a), 0])
        srg[f"{ang}deg_1m_de_T1"] = caja(q, 150)
    # perfiles a través de cada muro: media sobre la cara completa (|z| dentro de la altura interior) por capa de 10 cm
    # de profundidad, y sobre un parche de 1 × 1 m frente al blanco más cercano; error estándar entre corridas
    caras = {}
    zin = (zc > fz) & (zc < fz + ih_)
    srcs = [np.array(s["pos_mm"]) for s in geo["sources"]]
    for nombre, eje, signo in (("oeste", 0, -1), ("este", 0, 1), ("sur", 1, -1), ("norte", 1, 1)):
        semi = hx if eje == 0 else hy
        cen, tan = (xc, yc) if eje == 0 else (yc, xc)
        semi_t = hy if eje == 0 else hx
        # parche frente al blanco cuya proyección sobre la cara está más cerca
        pc = min(srcs, key=lambda q: abs(signo * semi - q[eje]))
        t0 = pc[1 - eje]
        cara_t = np.abs(tan) < semi_t
        parche_t = np.abs(tan - t0) <= 500
        parche_z = np.abs(zc - pc[2]) <= 500
        dd = np.arange(-300, wall + 1000, h)                    # cara inferior de la capa (mm); < 0 = sala
        filas = []
        for dmm in dd:
            i = int(np.argmin(np.abs(cen - signo * (semi + dmm + h / 2))))
            def media(mt, mz):
                if eje == 0:
                    return np.array([Hk[i][np.ix_(mt, mz)].mean() for Hk in H])
                return np.array([Hk[:, i, :][np.ix_(mt, mz)].mean() for Hk in H])
            vc, vp = media(cara_t, zin), media(parche_t, parche_z)
            se = lambda v: float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float("nan")
            filas.append(dict(prof_cm=(dmm + h / 2) / 10, cara=float(vc.mean() * k_uSvh), cara_se=se(vc) * k_uSvh,
                              parche=float(vp.mean() * k_uSvh), parche_se=se(vp) * k_uSvh))
        def ajuste(key):
            P = [(f["prof_cm"], f[key], f[key + "_se"]) for f in filas if 40 <= f["prof_cm"] <= wall / 10 - 10
                 and f[key] > 0 and f[key + "_se"] < 0.5 * f[key]]
            if len(P) < 4:
                return None, len(P)
            d_, v_, e_ = map(np.array, zip(*P))
            c = np.polyfit(d_, np.log(v_), 1, w=v_ / np.maximum(e_, 1e-30))
            return float(-c[0]), len(P)
        mu_c, n_c = ajuste("cara"); mu_p, n_p = ajuste("parche")
        ext = [f for f in filas if wall / 10 <= f["prof_cm"] <= wall / 10 + 30]
        interior = [f for f in filas if f["prof_cm"] < 0][-1]
        def comb(L, k):
            v = np.array([f[k] for f in L]); e = np.array([f[k + "_se"] for f in L])
            return float(v.mean()), float(np.sqrt(np.nansum(e ** 2)) / len(e))
        caras[nombre] = dict(
            blanco_frente=[round(float(t0)), round(float(pc[2]))],
            interior_cara_uSv_h=interior["cara"], interior_parche_uSv_h=interior["parche"],
            exterior_0_30cm_cara_uSv_h=comb(ext, "cara"), exterior_0_30cm_parche_uSv_h=comb(ext, "parche"),
            mu_cara_cm_inv=mu_c, tvl_cara_cm=(math.log(10) / mu_c if mu_c else None), capas_ajuste_cara=n_c,
            mu_parche_cm_inv=mu_p, tvl_parche_cm=(math.log(10) / mu_p if mu_p else None), capas_ajuste_parche=n_p,
            perfil=filas)
    integral_fuera_rel = float(Hm[fuera].sum() / Hm.sum())
    techo_ext = Hm[(np.abs(X) < hx) & (np.abs(Yg) < hy) & (Z > fz + ih_ + roof) & (Z < fz + ih_ + roof + 2 * h)]
    # activación de saturación
    mats = yaml.safe_load(open(ROOT / "materials.yaml"))
    regs = [leer_regiones(f.replace("_malla.bin", "_regiones.txt")) for f in mallas]
    # suma por corrida
    por_corrida = []
    for R in regs:
        acc = {}
        for k, v in R.items():
            kk = "hormigon" if isinstance(k, tuple) else k
            acc[kk] = acc.get(kk, 0) + v
        por_corrida.append(acc)
    act = {}
    for iso, x in XS.items():
        vals = []
        for acc in por_corrida:
            Rn = 0.0
            for reg in CANALES[iso]:
                if reg not in acc:
                    continue
                mname = "concrete_portland" if reg == "hormigon" else reg
                m = mats[mname]; w = m["wt"].get(x["el"], 0.0) / sum(m["wt"].values())
                Nd = m["density"] * w * NA * x["ab"] / x["A"]             # núcleos blanco por cm³
                th, ep = acc[reg][0], acc[reg][1]
                Rn += Nd * (x["s0"] * th + x["I0"] * ep / DU) * 1e-24    # reacciones por n fuente
            vals.append(Rn * S)
        vals = np.array(vals)
        act[iso] = dict(A_sat_Bq=float(vals.mean()), sigma_Bq=float(vals.std(ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else None,
                        canales=CANALES[iso], fuente_xs=x["fuente"], t12_d=x["t12_d"])
    flu = {}
    for acc in por_corrida:
        for k, v in acc.items():
            flu.setdefault(k, []).append(v)
    flu = {k: dict(termico_1v=float(np.mean([v[0] for v in L])), epitermico=float(np.mean([v[1] for v in L])),
                   rapido=float(np.mean([v[2] for v in L])), total=float(np.mean([v[3] for v in L]))) for k, L in flu.items()}
    out = dict(
        caso=str(d.relative_to(ROOT)), version_modelo=geo.get("version_modelo"), corridas=len(M), historias_totales=Ntot,
        fuente=dict(blancos=[s["name"] for s in geo["sources"]], I_uA_por_blanco=I_uA, rendimiento_n_por_p=Y, n_por_s_total=S,
                    puertas_HDPE=geo["puertas"]),
        bloque_hormigon=dict(interior_mm=[ix_, iy_, ih_], muro_mm=wall, losa_mm=slab, techo_mm=roof),
        srg_tabla8_uSv_h=srg, caras_plano_medio=caras, integral_dosis_fuera_bloque_rel=integral_fuera_rel,
        techo_exterior_medio_uSv_h=float(techo_ext.mean() * k_uSvh),
        activacion_saturacion=act, fluencia_por_region_cm_por_n=flu,
        malla=dict(h_mm=h, dims=[nx, ny, nz], origen_mm=[x0, y0, z0]))
    Mf = [m for m in M if "flux3g" in m]                            # corridas v1b con malla de fluencia
    if len(Mf) >= 2:
        Fl = np.array([m["flux3g"].astype(np.float64) / m["N"] / V for m in Mf])    # n/cm² por n fuente
        np.savez_compressed(d / "mallas_gemelo_f77.npz",
                            dose_h10_media=(Hm * k_uSvh).astype(np.float32), dose_h10_std=(Hs * k_uSvh).astype(np.float32),
                            flux3g_media=(Fl.mean(0) * S).astype(np.float32),
                            flux3g_std=(Fl.std(0, ddof=1) / math.sqrt(len(Mf)) * S).astype(np.float32),
                            corridas_dosis=len(M), corridas_flujo=len(Mf),
                            x_mm=xc, y_mm=yc, z_mm=zc, unidades=np.array(["µSv/h", "n/cm²/s (térmico <0.5 eV, epitérmico <0.1 MeV, rápido)"]))
    np.savez_compressed(d / "mapa_h10.npz", H_media_pSv_por_n=Hm.astype(np.float32), H_std_pSv_por_n=Hs.astype(np.float32),
                        Hn_media_pSv_por_n=Hn.mean(0).astype(np.float32), x_mm=xc, y_mm=yc, z_mm=zc, k_uSv_h_por_pSv_n=k_uSvh)
    json.dump(out, open(d / "resultado_mapa.json", "w"), indent=1, ensure_ascii=False, default=float)
    return out


if __name__ == "__main__":
    o = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "r*")
    print(json.dumps({k: o[k] for k in ("corridas", "historias_totales", "srg_tabla8_uSv_h", "integral_dosis_fuera_bloque_rel")}, default=float))
    for k, v in o["caras_plano_medio"].items():
        print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != "perfil"})
    print({k: f"{v['A_sat_Bq']:.3g} Bq" for k, v in o["activacion_saturacion"].items()})
