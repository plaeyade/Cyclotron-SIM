#!/usr/bin/env python3
"""Validación del GDML en Geant4 (geant4_pybind): solapes del navegador + transporte de geantinos.

No requiere datasets de física (geantino = solo transporte). Uso:
  echo n | python validate_g4.py out/mc/kiube_geant4.gdml [N_geantinos]
"""
import sys, json, math, collections, pathlib
import os, tempfile
if os.environ.get("KIUBE_G4_STUB", "1") == "1" and not os.environ.get("G4ENSDFSTATEDATA"):
    # modo geometría: Geant4 exige ENSDFSTATE aunque solo se transporten geantinos
    _d = tempfile.mkdtemp(); open(os.path.join(_d, "ENSDFSTATE.dat"), "w").close()
    os.environ["G4ENSDFSTATEDATA"] = _d
import geant4_pybind as g4

gdml = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
optics = json.load(open(pathlib.Path(gdml).parent / "sources.json"))

parser = g4.G4GDMLParser()
parser.SetOverlapCheck(False)
parser.Read(gdml, False)
world = parser.GetWorldVolume()
lv = world.GetLogicalVolume()

# ---------- 1) solapes (muestreo de puntos en superficie de cada volumen contra sus hermanos y la madre)
ov = []
for i in range(lv.GetNoDaughters()):
    pv = lv.GetDaughter(i)
    if pv.CheckOverlaps(2000, 0.001, False, 1):     # tolerancia 1 µm
        ov.append(pv.GetName())
print(f"[overlaps] {len(ov)} de {lv.GetNoDaughters()} volúmenes con solape > 1 µm: {ov[:20]}")


# ---------- 2) geantinos desde las cavidades de los blancos, isótropos
class Det(g4.G4VUserDetectorConstruction):
    def Construct(self):
        return world


class Phys(g4.G4VModularPhysicsList):
    """Solo transporte (G4VModularPhysicsList añade G4Transportation): no requiere datasets."""
    def ConstructParticle(self):
        g4.G4Geantino.GeantinoDefinition()

    def SetCuts(self):
        pass


class Gun(g4.G4VUserPrimaryGeneratorAction):
    def __init__(self):
        super().__init__()
        self.gun = g4.G4ParticleGun(1)
        self.gun.SetParticleDefinition(g4.G4Geantino.GeantinoDefinition())
        self.gun.SetParticleEnergy(10 * g4.MeV)
        self.k = 0

    def GeneratePrimaries(self, ev):
        s = optics["sources"][self.k % len(optics["sources"])]
        self.k += 1
        u = 2 * g4.G4UniformRand() - 1
        ph = 2 * math.pi * g4.G4UniformRand()
        st = math.sqrt(1 - u * u)
        self.gun.SetParticlePosition(g4.G4ThreeVector(*s["pos_mm"]) * g4.mm)
        self.gun.SetParticleMomentumDirection(g4.G4ThreeVector(st * math.cos(ph), st * math.sin(ph), u))
        self.gun.GeneratePrimaryVertex(ev)


path = collections.Counter()
first = collections.Counter()
stuck = [0]


class Step(g4.G4UserSteppingAction):
    def UserSteppingAction(self, step):
        pre = step.GetPreStepPoint()
        v = pre.GetTouchable().GetVolume()
        name = v.GetName() if v else "OUT"
        path[name] += step.GetStepLength() / g4.cm
        tr = step.GetTrack()
        if tr.GetCurrentStepNumber() == 1:
            first[name] += 1
        if tr.GetCurrentStepNumber() > 5000:
            stuck[0] += 1
            tr.SetTrackStatus(g4.G4TrackStatus.fStopAndKill)


rm = g4.G4RunManagerFactory.CreateRunManager(g4.G4RunManagerType.Serial)
rm.SetUserInitialization(Det())
rm.SetUserInitialization(Phys())
rm.SetUserAction(Gun())
rm.SetUserAction(Step())
rm.Initialize()
rm.BeamOn(N)
tot = sum(path.values())
print(f"[geantinos] {N} disparados; volumen de nacimiento: {dict(first)}")
print(f"[geantinos] tracks atascados (>5000 pasos): {stuck[0]}")
print("[geantinos] longitud de camino por volumen (top 15, cm):")
for k, v in path.most_common(15):
    print(f"   {k:28s} {v:12.1f}  ({100 * v / tot:5.2f} %)")
json.dump({"overlaps": ov, "n_daughters": lv.GetNoDaughters(), "stuck": stuck[0], "born_in": dict(first),
           "path_cm": dict(path)}, open(pathlib.Path(gdml).parent / "g4_validation.json", "w"), indent=1)
