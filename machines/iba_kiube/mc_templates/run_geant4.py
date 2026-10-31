#!/usr/bin/env python3
"""Geant4 (geant4_pybind ≥ 11.2) — neutrones de los blancos del Cyclone KIUBE sobre el GDML teselado.

  pip install geant4_pybind          # la 1.ª vez descarga los datasets (G4NDL etc., ~2 GB)
  python run_geant4.py --targets T1 T2 --events 100000 --physics QGSP_BIC_HP
  python run_geant4.py --geantino --events 2000      # prueba de geometría sin datasets

Salida: fluencia de neutrones por volumen (track-length / volumen), en 3 grupos de energía,
normalizada por neutrón fuente y por segundo de operación (sources.json).
"""
import argparse, json, math, bisect, pathlib, collections
import os, tempfile
if "--geantino" in __import__("sys").argv and not os.environ.get("G4ENSDFSTATEDATA"):
    # modo geometría: Geant4 exige ENSDFSTATE aunque solo se transporten geantinos
    _d = tempfile.mkdtemp(); open(os.path.join(_d, "ENSDFSTATE.dat"), "w").close()
    os.environ["G4ENSDFSTATEDATA"] = _d
import geant4_pybind as g4

HERE = pathlib.Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--targets", nargs="+", default=["T1", "T2"])
ap.add_argument("--events", type=int, default=100000)
ap.add_argument("--physics", default="QGSP_BIC_HP")
ap.add_argument("--geantino", action="store_true")
ap.add_argument("--threads", type=int, default=1)
a = ap.parse_args()

S = json.load(open(HERE / "sources.json"))
src = [s for s in S["sources"] if s["name"] in a.targets]
strength = sum(s["neutrons_per_s_if_active"] or 0 for s in src)
E = S["spectrum"]["E_MeV"]
w = S["spectrum"]["EdNdE_per_proton"]
o = sorted(range(len(E)), key=lambda i: E[i]); E = [E[i] for i in o]; w = [w[i] for i in o]
# CDF en log(E): dN = (E dN/dE) dlnE
cdf = [0.0]
for i in range(1, len(E)):
    cdf.append(cdf[-1] + 0.5 * (w[i] + w[i - 1]) * math.log(E[i] / E[i - 1]))
cdf = [c / cdf[-1] for c in cdf]


def sample_E():
    r = g4.G4UniformRand()
    i = max(1, bisect.bisect_left(cdf, r))
    f = (r - cdf[i - 1]) / max(cdf[i] - cdf[i - 1], 1e-30)
    return math.exp(math.log(E[i - 1]) + f * math.log(E[i] / E[i - 1]))


parser = g4.G4GDMLParser()
parser.Read(str(HERE / "kiube_geant4.gdml"), False)
WORLD = parser.GetWorldVolume()
VOL = {}
for i in range(WORLD.GetLogicalVolume().GetNoDaughters()):
    pv = WORLD.GetLogicalVolume().GetDaughter(i)
    VOL[pv.GetName()] = pv.GetLogicalVolume().GetSolid().GetCubicVolume() / g4.cm3


class Det(g4.G4VUserDetectorConstruction):
    def Construct(self):
        return WORLD


class GeantinoPhys(g4.G4VModularPhysicsList):
    """Solo transporte (G4VModularPhysicsList añade G4Transportation): no requiere datasets."""
    def ConstructParticle(self):
        g4.G4Geantino.GeantinoDefinition()

    def SetCuts(self):
        pass


class Gun(g4.G4VUserPrimaryGeneratorAction):
    def __init__(self):
        super().__init__()
        self.gun = g4.G4ParticleGun(1)
        part = g4.G4Geantino.GeantinoDefinition() if a.geantino else g4.G4Neutron.Definition()
        self.gun.SetParticleDefinition(part)
        self.k = 0

    def GeneratePrimaries(self, ev):
        s = src[self.k % len(src)]; self.k += 1
        u = 2 * g4.G4UniformRand() - 1; ph = 2 * math.pi * g4.G4UniformRand(); st = math.sqrt(1 - u * u)
        self.gun.SetParticleEnergy(sample_E() * g4.MeV)
        self.gun.SetParticlePosition(g4.G4ThreeVector(*s["pos_mm"]) * g4.mm)
        self.gun.SetParticleMomentumDirection(g4.G4ThreeVector(st * math.cos(ph), st * math.sin(ph), u))
        self.gun.GeneratePrimaryVertex(ev)


GROUPS = (0.625e-6, 0.1)                 # MeV: térmico < 0.625 eV < epitérmico < 100 keV < rápido
L = collections.defaultdict(lambda: [0.0, 0.0, 0.0])
target_pdg = 0 if a.geantino else 2112


class Step(g4.G4UserSteppingAction):
    def UserSteppingAction(self, step):
        tr = step.GetTrack()
        if tr.GetDefinition().GetPDGEncoding() != target_pdg:
            return
        pre = step.GetPreStepPoint()
        v = pre.GetTouchable().GetVolume()
        if v is None:
            return
        e = pre.GetKineticEnergy() / g4.MeV
        gi = 0 if e < GROUPS[0] else (1 if e < GROUPS[1] else 2)
        L[v.GetName()][gi] += step.GetStepLength() / g4.cm


rm = g4.G4RunManagerFactory.CreateRunManager(
    g4.G4RunManagerType.Serial if a.threads == 1 else g4.G4RunManagerType.MT)
if a.threads > 1:
    rm.SetNumberOfThreads(a.threads)
rm.SetUserInitialization(Det())
if a.geantino:
    rm.SetUserInitialization(GeantinoPhys())
else:
    rm.SetUserInitialization(g4.G4PhysListFactory().GetReferencePhysList(a.physics))
rm.SetUserAction(Gun()); rm.SetUserAction(Step())
rm.Initialize()
rm.BeamOn(a.events)

res = {"events": a.events, "targets": a.targets, "neutrons_per_s": strength, "physics": "geantino" if a.geantino else a.physics,
       "fluence_per_source_n_cm2": {}, "fluence_rate_n_cm2_s": {}}
for name, l3 in sorted(L.items(), key=lambda kv: -sum(kv[1])):
    vol = VOL.get(name)
    if not vol:
        continue
    f = [x / vol / a.events for x in l3]
    res["fluence_per_source_n_cm2"][name] = {"thermal": f[0], "epithermal": f[1], "fast": f[2]}
    res["fluence_rate_n_cm2_s"][name] = sum(f) * strength
json.dump(res, open("kiube_geant4_results.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ("events", "targets", "physics")}))
for k, v in list(res["fluence_rate_n_cm2_s"].items())[:12]:
    print(f"  {k:28s} {v:.3e} n/cm²/s")
