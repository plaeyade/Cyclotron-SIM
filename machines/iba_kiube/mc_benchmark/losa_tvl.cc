// losa_tvl — transmisión de H*(10) (n + γ, ICRP 74) a través de una losa semi-infinita de material de blindaje para el
// espectro 18O(p,n) BON21 de kiube (haz delgado normal; la integral sobre el plano a profundidad z equivale, por
// invariancia de traslación, a la H*(10) en un punto bajo un haz ancho uniforme de fluencia unitaria).
// Tally: longitud de traza ponderada con h(E) por capa de 10 cm DENTRO de la losa (pSv·cm² por n incidente).
// Reducción de varianza: importancia geométrica por capa (f^k) en mundo paralelo, n y γ.
// Uso: losa_tvl material.txt espectro_config.txt N semilla espesor_cm f salida.txt
//   material.txt: "rho <g/cm3>" y "wt El frac El frac ..."
#include "G4RunManagerFactory.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserStackingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4ParticleGun.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

static const double nE[] = {1e-9, 1e-8, 2.53e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4,
                            5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9,
                            1.0, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 15.0, 16.0, 18.0, 20.0};
static const double nH[] = {6.60, 9.00, 10.6, 12.9, 13.5, 13.6, 13.3, 12.9, 12.0, 11.3, 10.6, 9.90, 9.40, 8.90, 8.30,
                            7.90, 7.70, 8.00, 10.5, 16.6, 23.7, 41.1, 60.0, 88.0, 132, 170, 233, 322, 375, 400, 416, 425,
                            420, 412, 408, 405, 400, 405, 409, 420, 440, 480, 520, 540, 555, 570, 600};
static const double gE[] = {0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
                            1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
static const double gH[] = {0.061, 0.83, 1.05, 0.81, 0.64, 0.55, 0.51, 0.53, 0.61, 0.89, 1.20, 1.80, 2.38, 2.93, 3.44,
                            4.38, 5.20, 6.90, 8.60, 11.1, 13.4, 15.5, 17.6, 21.6, 25.6};
static double loglog(const double* x, const double* y, int n, double e) {
  if (e <= x[0]) return (n == 47) ? y[0] : y[0] * e / x[0];
  if (e >= x[n - 1]) return y[n - 1];
  int i = 1; while (x[i] < e) ++i;
  double t = std::log(e / x[i - 1]) / std::log(x[i] / x[i - 1]);
  return std::exp(std::log(y[i - 1]) + t * std::log(y[i] / y[i - 1]));
}

static double RHO = 2.35, TH = 200.0, F = 2.0; static const double DZ = 10.0;   // cm
static std::vector<std::pair<std::string, double>> WT;
static std::vector<double> SE, CDF;
static int NL;
static std::vector<double> hN, hG, hN2, hG2, thF;          // por capa: suma, suma de cuadrados por historia, fluencia térmica
static std::vector<double> evN, evG;                        // acumuladores del evento
static G4VPhysicalVolume* gSlab = nullptr;

class Geo : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* mat = new G4Material("blindaje", RHO * g / cm3, (G4int)WT.size());
    double tot = 0; for (auto& e : WT) tot += e.second;
    for (auto& e : WT) mat->AddElement(nist->FindOrBuildElement(e.first), e.second / tot);
    auto* vac = nist->FindOrBuildMaterial("G4_Galactic");
    const double L = 20 * m, Tz = TH * cm;
    auto* wl = new G4LogicalVolume(new G4Box("w", L, L, Tz + 1 * m), vac, "w");
    auto* W = new G4PVPlacement(nullptr, {}, wl, "w", nullptr, false, 0);
    auto* sl = new G4LogicalVolume(new G4Box("s", L - 1 * mm, L - 1 * mm, Tz / 2), mat, "s");
    gSlab = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, Tz / 2), sl, "s", wl, false, 0);
    return W;
  }
};
class Imp : public G4VUserParallelWorld {
 public:
  Imp() : G4VUserParallelWorld("imp") {}
  void Construct() override {
    auto* gl = GetWorld()->GetLogicalVolume();
    for (int k = 0; k < NL; ++k) {
      auto* L = new G4LogicalVolume(new G4Box("c", 20 * m - 2 * mm, 20 * m - 2 * mm, DZ * cm / 2), nullptr, "c");
      cells.push_back(new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (k + 0.5) * DZ * cm), L, "c" + std::to_string(k), gl, false, k));
    }
    w = GetWorld();
  }
  void CreateStore() {
    auto* st = G4IStore::GetInstance(GetName());
    st->AddImportanceGeometryCell(1.0, *w, 0);
    for (int k = 0; k < NL; ++k) st->AddImportanceGeometryCell(std::pow(F, k), *cells[k], k);
  }
  std::vector<G4VPhysicalVolume*> cells; G4VPhysicalVolume* w = nullptr;
};
class Gen : public G4VUserPrimaryGeneratorAction {
 public:
  Gen() { g = new G4ParticleGun(1); g->SetParticleDefinition(G4Neutron::Definition()); }
  void GeneratePrimaries(G4Event* ev) override {
    double r = G4UniformRand(); size_t i = 1; while (i < CDF.size() - 1 && CDF[i] < r) ++i;
    double f = (r - CDF[i - 1]) / std::max(CDF[i] - CDF[i - 1], 1e-30);
    g->SetParticleEnergy(std::exp(std::log(SE[i - 1]) + f * std::log(SE[i] / SE[i - 1])) * MeV);
    g->SetParticlePosition(G4ThreeVector(0, 0, -1 * mm)); g->SetParticleMomentumDirection({0, 0, 1});
    g->GeneratePrimaryVertex(ev);
  }
  G4ParticleGun* g;
};
class Step : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* st) override {
    const G4Track* tr = st->GetTrack(); auto* d = tr->GetDefinition();
    bool n = d == G4Neutron::Definition(), gm = d == G4Gamma::Definition();
    if (!n && !gm) return;
    auto* a = st->GetPreStepPoint(); if (a->GetPhysicalVolume() != gSlab) return;
    double L = st->GetStepLength() / cm; if (L <= 0) return;
    double z0 = a->GetPosition().z() / cm, z1 = st->GetPostStepPoint()->GetPosition().z() / cm;
    double E = a->GetKineticEnergy() / MeV, w = a->GetWeight();
    double h = n ? loglog(nE, nH, 47, E) : loglog(gE, gH, 25, E);
    // reparte la traza entre capas según z
    double lo = std::min(z0, z1), hi = std::max(z0, z1), dz = hi - lo;
    int k0 = std::max(0, (int)(lo / DZ)), k1 = std::min(NL - 1, (int)(hi / DZ));
    for (int k = k0; k <= k1; ++k) {
      double frac = dz > 1e-9 ? (std::min(hi, (k + 1) * DZ) - std::max(lo, k * DZ)) / dz : 1.0;
      if (frac <= 0) continue;
      double v = w * h * L * frac / DZ;
      if (n) { evN[k] += v; if (E < 5e-7) thF[k] += w * L * frac / DZ; } else evG[k] += v;
    }
  }
};
class Stack : public G4UserStackingAction {
 public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
    auto* p = t->GetDefinition();
    return (p == G4Neutron::Definition() || p == G4Gamma::Definition()) ? fUrgent : fKill;
  }
  void NewStage() override {}
  void PrepareNewEvent() override {
    // cierre del evento anterior: acumula sumas y cuadrados (error por historia)
    for (int k = 0; k < NL; ++k) { hN[k] += evN[k]; hN2[k] += evN[k] * evN[k]; hG[k] += evG[k]; hG2[k] += evG[k] * evG[k]; }
    std::fill(evN.begin(), evN.end(), 0.0); std::fill(evG.begin(), evG.end(), 0.0);
  }
};
class Act : public G4VUserActionInitialization {
 public:
  void Build() const override { SetUserAction(new Gen); SetUserAction(new Step); SetUserAction(new Stack); }
};

int main(int argc, char** argv) {
  if (argc < 8) { std::cerr << "uso: losa_tvl material.txt espectro_config.txt N semilla espesor_cm f salida\n"; return 1; }
  { std::ifstream in(argv[1]); std::string l;
    while (std::getline(in, l)) { std::istringstream s(l); std::string k; s >> k;
      if (k == "rho") s >> RHO; else if (k == "wt") { std::string e; double v; while (s >> e >> v) WT.push_back({e, v}); } } }
  { std::ifstream in(argv[2]); std::string l; std::vector<double> w;
    while (std::getline(in, l)) { std::istringstream s(l); std::string k; s >> k;
      if (k == "spectrum") { int n; s >> n; for (int i = 0; i < n; ++i) { std::getline(in, l); std::istringstream t(l); double e, y; t >> e >> y; SE.push_back(e); w.push_back(y); } } }
    CDF.assign(SE.size(), 0.0);
    for (size_t i = 1; i < SE.size(); ++i) CDF[i] = CDF[i - 1] + 0.5 * (w[i] + w[i - 1]) * std::log(SE[i] / SE[i - 1]);
    for (auto& c : CDF) c /= CDF.back(); }
  long N = std::atol(argv[3]); G4Random::setTheSeed(std::atol(argv[4])); TH = std::atof(argv[5]); F = std::atof(argv[6]);
  NL = (int)std::lround(TH / DZ);
  hN.assign(NL, 0); hG.assign(NL, 0); hN2.assign(NL, 0); hG2.assign(NL, 0); thF.assign(NL, 0); evN.assign(NL, 0); evG.assign(NL, 0);
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* geo = new Geo; auto* imp = new Imp; geo->RegisterParallelWorld(imp);
  rm->SetUserInitialization(geo);
  G4PhysListFactory fac; auto* pl = fac.GetReferencePhysList("QGSP_BIC_HP"); pl->SetVerboseLevel(0);
  auto* sN = new G4GeometrySampler("imp", "neutron"); sN->SetParallel(true);
  auto* sG = new G4GeometrySampler("imp", "gamma"); sG->SetParallel(true);
  pl->RegisterPhysics(new G4ImportanceBiasing(sN, "imp")); pl->RegisterPhysics(new G4ImportanceBiasing(sG, "imp"));
  pl->RegisterPhysics(new G4ParallelWorldPhysics("imp"));
  rm->SetUserInitialization(pl); rm->SetUserInitialization(new Act);
  rm->Initialize(); imp->CreateStore();
  rm->BeamOn((G4int)N);
  for (int k = 0; k < NL; ++k) { hN[k] += evN[k]; hN2[k] += evN[k] * evN[k]; hG[k] += evG[k]; hG2[k] += evG[k] * evG[k]; }
  std::ofstream o(argv[7]);
  o << "# losa_tvl rho " << RHO << " N " << N << " f " << F << " seed " << argv[4] << "\n# z_cm Hn Hn_se Hg Hg_se phi_term  [pSv cm2 / n incidente; cm2 / n]\n";
  for (int k = 0; k < NL; ++k) {
    double mn = hN[k] / N, mg = hG[k] / N;
    double sn = std::sqrt(std::max(0.0, hN2[k] / N - mn * mn) / N), sg = std::sqrt(std::max(0.0, hG2[k] / N - mg * mg) / N);
    o << std::setprecision(6) << (k + 0.5) * DZ << " " << mn << " " << sn << " " << mg << " " << sg << " " << thF[k] / N << "\n";
  }
  delete rm; return 0;
}
