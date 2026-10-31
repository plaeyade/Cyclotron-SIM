"""Monte Carlo de partículas de prueba (TPMC) en flujo molecular — la misma física que Molflow+.

Moléculas independientes (régimen molecular, Kn ≫ 1), emisión y re-emisión difusa con ley del
coseno, probabilidad de captura (sticking) por faceta para las bombas, fuentes de gas por faceta
(desgasificación específica en mbar·l/s/cm² o caudal fijo). Trazado de rayos con Embree
(trimesh + embreex) sobre el STL de superficies del vacío (kiube_vacio_cm.stl, LOD 'vacuum').

Presión en una faceta a partir de la tasa de impacto (equilibrio isotrópico, Maxwell):
    J = n·v̄/4  →  P = n·k·T = 4·J·k·T / v̄
Con Q [mbar·l/s] repartido en N partículas de prueba, cada impacto representa Q/N; k·T se cancela:
    P [mbar] = 4 · hits · (Q·1e-3 m³·100 Pa / N) / (A [m²] · v̄ [m/s]) · 1e-2
El flujo molecular es lineal: la respuesta a varias fuentes es la suma de las respuestas, y la
evolución lenta (horas) del pump-down es cuasi-estática (V/S ~ 0.1 s ≪ t) → P(t) = Σ_i q_i(t)·R_i.

Referencias de la física: Kersevan & Pons, JVST A 27 (2009) 1017 (Molflow+); Clausing (1932);
Chiggiato, CAS 2017, arXiv:2006.07124 (desgasificación).
"""
from __future__ import annotations

import json, math, pathlib
import numpy as np

K_B = 1.380649e-23
AMU = 1.66053907e-27


def vbar(mass_g_mol, T=293.15):
    """Velocidad media de Maxwell [m/s]."""
    return math.sqrt(8 * K_B * T / (math.pi * mass_g_mol * AMU))


def sticking_for_speed(S_l_s, area_cm2, mass_g_mol, T=293.15):
    """Probabilidad de captura que da la velocidad de bombeo S sobre el área A: S = s·A·v̄/4."""
    return min(1.0, 4 * S_l_s * 1e-3 / (area_cm2 * 1e-4 * vbar(mass_g_mol, T)))


def _cosine_dirs(normals, rng):
    """Direcciones con ley del coseno alrededor de cada normal (vectorizado)."""
    n = len(normals)
    u1, u2 = rng.random(n), rng.random(n)
    st, ct = np.sqrt(u1), np.sqrt(1 - u1)               # sin θ, cos θ con pdf ∝ cos θ sin θ
    ph = 2 * np.pi * u2
    # base ortonormal por normal
    a = np.where(np.abs(normals[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    t1 = np.cross(normals, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(normals, t1)
    return (st * np.cos(ph))[:, None] * t1 + (st * np.sin(ph))[:, None] * t2 + ct[:, None] * normals


def _sample_on_tris(V, F, tri_idx, rng):
    a, b, c = V[F[tri_idx, 0]], V[F[tri_idx, 1]], V[F[tri_idx, 2]]
    r1, r2 = rng.random(len(tri_idx)), rng.random(len(tri_idx))
    s = np.sqrt(r1)
    return (1 - s)[:, None] * a + (s * (1 - r2))[:, None] * b + (s * r2)[:, None] * c


class TPMC:
    def __init__(self, mesh, sticking, seed=1):
        """mesh: trimesh.Trimesh (unidades arbitrarias, consistentes; aquí cm). sticking: array por faceta."""
        import trimesh  # noqa: F401
        self.m = mesh
        self.V = np.asarray(mesh.vertices, float)
        self.F = np.asarray(mesh.faces)
        self.N = np.asarray(mesh.face_normals, float)
        self.A = np.asarray(mesh.area_faces, float)
        self.s = np.asarray(sticking, float)
        self.rng = np.random.default_rng(seed)
        self.eps = 1e-4 * float(np.sqrt(np.median(self.A)))
        self.ray = mesh.ray                                 # Embree si embreex está instalado
        self.uses_embree = "embree" in type(self.ray).__module__

    def run(self, source_weights, n_particles=20000, batch=20000, max_bounces=200000):
        """Lanza n_particles desde facetas fuente (pesos ∝ caudal). Devuelve impactos y capturas por faceta."""
        w = np.asarray(source_weights, float)
        src = np.flatnonzero(w > 0)
        p = w[src] / w[src].sum()
        hits = np.zeros(len(self.F)); absorbed = np.zeros(len(self.F))
        lost = 0; bounces_total = 0
        done = 0
        while done < n_particles:
            nb = min(batch, n_particles - done)
            tri = self.rng.choice(src, size=nb, p=p)
            pos = _sample_on_tris(self.V, self.F, tri, self.rng)
            nrm = self.N[tri] * np.sign(self._inward(tri))[:, None]
            d = _cosine_dirs(nrm, self.rng)
            pos = pos + self.eps * nrm
            alive = np.ones(nb, bool)
            for _ in range(max_bounces):
                idx = np.flatnonzero(alive)
                if not len(idx):
                    break
                loc, iray, itri = self.ray.intersects_location(pos[idx], d[idx], multiple_hits=False)
                miss = np.ones(len(idx), bool); miss[iray] = False
                if miss.any():                               # escapes por aberturas residuales
                    lost += int(miss.sum()); alive[idx[miss]] = False
                pid = idx[iray]
                np.add.at(hits, itri, 1.0)
                bounces_total += len(itri)
                stick = self.rng.random(len(itri)) < self.s[itri]
                np.add.at(absorbed, itri[stick], 1.0)
                alive[pid[stick]] = False
                keep = ~stick
                pid, itri, loc = pid[keep], itri[keep], loc[keep]
                nn = self.N[itri]
                flip = np.einsum("ij,ij->i", nn, d[pid]) > 0     # normal hacia el lado de llegada
                nn = np.where(flip[:, None], -nn, nn)
                d[pid] = _cosine_dirs(nn, self.rng)
                pos[pid] = loc + self.eps * nn
            lost += int(alive.sum())
            done += nb
        return dict(hits=hits, absorbed=absorbed, n=n_particles, lost=lost,
                    mean_bounces=bounces_total / n_particles)

    def run_tof(self, source_weights, detector, t_edges, mass_g_mol, T=293.15, n_particles=20000,
                batch=20000, max_bounces=200000, length_scale_m=1e-2):
        """Modo dependiente del tiempo (como Molflow 'time-dependent'): cada partícula de prueba sale en
        t = 0 y lleva su tiempo de vuelo; la rapidez se sorteа en cada (re)emisión de la Maxwelliana de
        flujo (∝ v³ e^(−v²/a²), acomodación térmica completa, como Molflow). Devuelve el histograma de
        impactos sobre las facetas 'detector' por intervalo de tiempo (respuesta al impulso)."""
        a = math.sqrt(2 * 1.380649e-23 * T / (mass_g_mol * 1.66053907e-27))   # velocidad más probable
        w = np.asarray(source_weights, float); src = np.flatnonzero(w > 0); p = w[src] / w[src].sum()
        det = np.asarray(detector, bool)
        H = np.zeros(len(t_edges) - 1); tmax = t_edges[-1]

        def speed(n):                                         # v = a·sqrt(G), G ~ Gamma(2, 1)
            return a * np.sqrt(self.rng.gamma(2.0, 1.0, n))
        done = 0
        while done < n_particles:
            nb = min(batch, n_particles - done)
            tri = self.rng.choice(src, size=nb, p=p)
            pos = _sample_on_tris(self.V, self.F, tri, self.rng)
            nrm = self.N[tri]
            d = _cosine_dirs(nrm, self.rng); pos = pos + self.eps * nrm
            v = speed(nb); t = np.zeros(nb)
            alive = np.ones(nb, bool)
            for _ in range(max_bounces):
                idx = np.flatnonzero(alive)
                if not len(idx):
                    break
                loc, iray, itri = self.ray.intersects_location(pos[idx], d[idx], multiple_hits=False)
                miss = np.ones(len(idx), bool); miss[iray] = False
                alive[idx[miss]] = False
                pid = idx[iray]
                t[pid] += np.linalg.norm(loc - pos[pid], axis=1) * length_scale_m / v[pid]
                late = t[pid] > tmax
                alive[pid[late]] = False
                pid, itri, loc = pid[~late], itri[~late], loc[~late]
                m = det[itri]
                if m.any():
                    H += np.histogram(t[pid[m]], bins=t_edges)[0]
                stick = self.rng.random(len(itri)) < self.s[itri]
                alive[pid[stick]] = False
                keep = ~stick
                pid, itri, loc = pid[keep], itri[keep], loc[keep]
                nn = self.N[itri]
                flip = np.einsum("ij,ij->i", nn, d[pid]) > 0
                nn = np.where(flip[:, None], -nn, nn)
                d[pid] = _cosine_dirs(nn, self.rng); pos[pid] = loc + self.eps * nn
                v[pid] = speed(len(pid))
            done += nb
        return dict(H=H, n=n_particles, t_edges=np.asarray(t_edges))

    def _inward(self, tri):
        # las normales del STL de la cámara ya apuntan al gas (vacuum.classify_vacuum_surfaces);
        # para mallas de prueba se usa la misma convención
        return np.ones(len(tri))

    def pressure(self, res, Q_mbar_l_s, mass_g_mol, T=293.15, area_scale_m2=1e-4):
        """Presión por faceta [mbar] para un caudal total Q de esa corrida."""
        vb = vbar(mass_g_mol, T)
        A_m2 = self.A * area_scale_m2
        per = Q_mbar_l_s * 1e-3 * 100.0 / res["n"]          # Pa·m³/s por partícula de prueba
        return 4 * res["hits"] * per / (A_m2 * vb) * 1e-2


def load_gas_mesh(path):
    """Carga el STL del vacío, fusiona vértices y orienta todas las normales HACIA EL GAS sin
    reordenar facetas (los rangos de índices por grupo siguen válidos):
    cada cáscara cerrada se hace consistente (fix_normals → hacia afuera de su volumen); la cáscara
    que encierra el gas (mayor volumen) se invierte; las interiores (láminas) quedan hacia afuera."""
    import trimesh
    m0 = trimesh.load(path, process=False)
    m = trimesh.Trimesh(m0.vertices, m0.faces, process=True, validate=False)
    assert len(m.faces) == len(m0.faces)
    labels = trimesh.graph.connected_component_labels(m.face_adjacency, node_count=len(m.faces))
    F = m.faces.copy()
    vols = {}
    for lab in np.unique(labels):
        idx = np.flatnonzero(labels == lab)
        sub = trimesh.Trimesh(m.vertices, m.faces[idx], process=False)
        trimesh.repair.fix_winding(sub); trimesh.repair.fix_inversion(sub)
        F[idx] = sub.faces
        vols[lab] = abs(sub.volume) if sub.is_watertight else 0.0
    outer = max(vols, key=vols.get)
    idx = np.flatnonzero(labels == outer)
    F[idx] = F[idx][:, ::-1]                                    # cavidad: normales hacia dentro (gas)
    out = trimesh.Trimesh(m.vertices, F, process=False)
    out.metadata["shells"] = len(vols)
    out.metadata["watertight"] = bool(m.is_watertight)
    return out


# ------------------------------------------------------------------ mallas de prueba (benchmarks)
def tube_mesh(R=1.0, L=2.0, n=96):
    """Tubo abierto con tapas = facetas de entrada/salida (normales hacia dentro). cm."""
    import trimesh
    cyl = trimesh.creation.cylinder(radius=R, height=L, sections=n)
    cyl.apply_translation([0, 0, L / 2])
    cyl.invert()                                          # normales hacia el interior (gas)
    z = cyl.triangles_center[:, 2]
    entrada = np.isclose(z, 0.0, atol=1e-9) & (np.abs(cyl.face_normals[:, 2]) > 0.99)
    salida = np.isclose(z, L, atol=1e-9) & (np.abs(cyl.face_normals[:, 2]) > 0.99)
    return cyl, entrada, salida


def box_with_pump(side_cm, pump_r_cm, n=64):
    """Cubo cerrado de lado L con un disco de bombeo centrado en una cara. Normales hacia dentro."""
    import trimesh
    b = trimesh.creation.box(extents=[side_cm] * 3)
    b = b.subdivide().subdivide().subdivide()
    b.invert()
    c = b.triangles_center
    pump = (np.isclose(c[:, 2], side_cm / 2)) & (np.hypot(c[:, 0], c[:, 1]) < pump_r_cm)
    return b, pump


# Clausing (1932) / Cole (1976): probabilidad de transmisión de un tubo cilíndrico, L/R
CLAUSING = {1.0: 0.6720, 2.0: 0.5136, 10.0: 0.1910}
