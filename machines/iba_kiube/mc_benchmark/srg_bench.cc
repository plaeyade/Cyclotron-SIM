// srg_bench — benchmark de blindaje del Cyclone KIUBE frente al SRG de IBA (MID 55672 Rev. G, Tablas 8 y 9).
//
// Geometría: ciclotrón teselado (GDML del LOD neutrónico, conforme) dentro del búnker estándar del SRG
// (interior 4000 × 4000 × 2500 mm, hormigón ρ = 2.35 g/cm³, muros y techo de espesor 't').
// Fuente: neutrones de 18O(p,n) desde la cavidad de los blancos F18 activos (doble haz: T1 y T5),
//         espectro BON21 (FLUKA, Bonvin et al. 2021), emisión isótropa. Las gammas de captura/inelásticas
//         se generan en el transporte (QGSP_BIC_HP, datos G4NDL 4.7.1).
// Estimador: longitud de traza en volúmenes de aire × coeficiente H*(10)/Φ de ICRP 74 (n y γ), por partícula fuente.
// Reducción de varianza: importancia geométrica (división/ruleta) en capas de 10 cm del muro frontal (+y),
//         mundo paralelo, factor 'f' por capa, para neutrones y gammas.
//
// Uso: srg_bench config.txt kiube.gdml t_muro_cm n_eventos semilla salida.txt [f=2.0]
#include "G4RunManagerFactory.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4GeometryCell.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4GDMLParser.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4ParticleGun.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>

// ---------------------------------------------------------------- coeficientes ICRP 74 H*(10)/Φ [pSv cm²]
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
  if (e <= x[0]) return (n == 47) ? y[0] : y[0] * e / x[0];   // n: térmico constante; γ < 10 keV ∝ E (despreciable)
  if (e >= x[n - 1]) return y[n - 1];
  int i = 1; while (x[i] < e) ++i;
  double t = std::log(e / x[i - 1]) / std::log(x[i] / x[i - 1]);
  return std::exp(std::log(y[i - 1]) + t * std::log(y[i] / y[i - 1]));
}

// ---------------------------------------------------------------- configuración
struct Src { std::string name; G4ThreeVector pos, dir; };
struct Cfg {
  double floor_z = -1120, ix = 4000, iy = 4000, ih = 2500, rho = 2.35, t = 2300, f = 2.0;
  std::vector<std::pair<std::string, double>> wt;
  std::vector<Src> src;
  std::vector<double> E, cdf;
};
static Cfg C;
static void readCfg(const char* fn) {
  std::ifstream in(fn); std::string line;
  std::vector<double> w;
  while (std::getline(in, line)) {
    std::istringstream s(line); std::string k; s >> k;
    if (k == "floor_z") s >> C.floor_z;
    else if (k == "vault_inner") s >> C.ix >> C.iy >> C.ih;
    else if (k == "concrete_density") s >> C.rho;
    else if (k == "concrete_wt") { std::string e; double v; while (s >> e >> v) C.wt.push_back({e, v}); }
    else if (k == "source") { Src x; double a, b, c, d, e; s >> x.name >> a >> b >> c >> d >> e;
      x.pos = G4ThreeVector(a, b, c); x.dir = G4ThreeVector(d, e, 0); C.src.push_back(x); }
    else if (k == "spectrum") {
      int n; s >> n;
      for (int i = 0; i < n; ++i) { std::getline(in, line); std::istringstream t(line); double e, y; t >> e >> y; C.E.push_back(e); w.push_back(y); }
    }
  }
  C.cdf.assign(C.E.size(), 0.0);                       // CDF en ln E de (E dN/dE)
  for (size_t i = 1; i < C.E.size(); ++i) C.cdf[i] = C.cdf[i - 1] + 0.5 * (w[i] + w[i - 1]) * std::log(C.E[i] / C.E[i - 1]);
  for (auto& c : C.cdf) c /= C.cdf.back();
}

// ---------------------------------------------------------------- detectores (volúmenes de aire)
struct Det { std::string name; G4ThreeVector c, half; G4VPhysicalVolume* pv = nullptr; double vol_cm3 = 0; };
static std::vector<Det> D;
static std::map<const G4VPhysicalVolume*, int> DIDX;
static const G4VPhysicalVolume* gVault = nullptr;
static double gLayN[40] = {0}, gLayW[40] = {0};

class Geo : public G4VUserDetectorConstruction {
 public:
  explicit Geo(const char* gdml) : fG(gdml) {}
  G4VPhysicalVolume* Construct() override {
    G4GDMLParser p; p.SetOverlapCheck(false); p.Read(fG, false);
    G4VPhysicalVolume* W = p.GetWorldVolume();
    G4LogicalVolume* WL = W->GetLogicalVolume();
    G4Material* air = WL->GetMaterial();
    auto* nist = G4NistManager::Instance();
    auto* conc = new G4Material("concrete_srg", C.rho * g / cm3, (G4int)C.wt.size());
    double tot = 0; for (auto& e : C.wt) tot += e.second;
    for (auto& e : C.wt) conc->AddElement(nist->FindOrBuildElement(e.first), e.second / tot);
    const double t = C.t, hx = C.ix / 2, hy = C.iy / 2, zf = C.floor_z, fl = 1000.0;
    auto* out = new G4Box("vo", hx + t, hy + t, (C.ih + t + fl) / 2);          // muros 't', techo 't', losa 1 m
    auto* inn = new G4Box("vi", hx, hy, C.ih / 2);
    G4ThreeVector co(0, 0, zf - fl + (C.ih + t + fl) / 2), ci(0, 0, zf + C.ih / 2);
    auto* shell = new G4SubtractionSolid("vault", out, inn, nullptr, ci - co);
    gVault = new G4PVPlacement(nullptr, co, new G4LogicalVolume(shell, conc, "vault"), "vault", WL, false, 0, true);
    const Src& s = C.src[0];
    double xh = s.pos.x() + (hy - s.pos.y()) * s.dir.x() / s.dir.y();          // eje del haz de T1 en y = +hy
    D = {{"frente_0deg_interior", {xh, hy - 150, 0}, {150, 150, 150}},
         {"lateral_90deg_interior", {hx - 150, s.pos.y(), 0}, {150, 150, 150}},
         {"opuesto_180deg_interior", {-xh, -hy + 150, 0}, {150, 150, 150}},
         {"a_1m_del_blanco_0deg", s.pos + 1000 * s.dir, {150, 150, 150}},
         {"exterior_tras_muro", {xh, hy + t + 100, 0}, {500, 100, 500}}};
    for (size_t i = 0; i < D.size(); ++i) {
      auto& d = D[i];
      auto* b = new G4Box(d.name, d.half.x(), d.half.y(), d.half.z());
      d.pv = new G4PVPlacement(nullptr, d.c, new G4LogicalVolume(b, air, d.name), d.name, WL, false, 0, true);
      d.vol_cm3 = 8 * d.half.x() * d.half.y() * d.half.z() / 1000.0;
      DIDX[d.pv] = (int)i;
    }
    {   // resto de la cara exterior del muro frontal (cara completa menos el detector central): promedio de cara
      const Det& dc = D.back();
      G4ThreeVector cb(0, hy + t + 100, co.z());
      auto* big = new G4Box("cara", hx + t, 100, (C.ih + t + fl) / 2);
      auto* hole = new G4Box("hueco", dc.half.x(), dc.half.y() + 1, dc.half.z());
      auto* frame = new G4SubtractionSolid("exterior_resto_cara", big, hole, nullptr, dc.c - cb);
      Det r; r.name = "exterior_resto_cara"; r.c = cb; r.half = G4ThreeVector(hx + t, 100, (C.ih + t + fl) / 2);
      r.pv = new G4PVPlacement(nullptr, cb, new G4LogicalVolume(frame, air, r.name), r.name, WL, false, 0, true);
      r.vol_cm3 = (8 * r.half.x() * r.half.y() * r.half.z() - 8 * dc.half.x() * dc.half.y() * dc.half.z()) / 1000.0;
      D.push_back(r); DIDX[r.pv] = (int)D.size() - 1;
    }
    return W;
  }
 private:
  const char* fG;
};

// mundo paralelo de importancias: capas de 10 cm en el muro frontal (+y) y una capa final hasta el borde
class Imp : public G4VUserParallelWorld {
 public:
  Imp() : G4VUserParallelWorld("imp") {}
  void Construct() override {
    G4VPhysicalVolume* gw = GetWorld();
    G4LogicalVolume* gl = gw->GetLogicalVolume();
    auto* wb = (G4Box*)gl->GetSolid();
    double X = wb->GetXHalfLength(), Y = wb->GetYHalfLength(), Z = wb->GetZHalfLength();
    double y0 = C.iy / 2, step = 100.0;
    int n = (int)std::ceil(C.t / step - 1e-9);
    for (int k = 0; k <= n; ++k) {
      double a = y0 + k * step, b = (k < n) ? std::min(y0 + (k + 1) * step, y0 + C.t) : Y;
      auto* box = new G4Box("cap" + std::to_string(k), X, (b - a) / 2, Z);
      auto* pv = new G4PVPlacement(nullptr, G4ThreeVector(0, (a + b) / 2, 0),
                                   new G4LogicalVolume(box, nullptr, "cap" + std::to_string(k)),
                                   "cap" + std::to_string(k), gl, false, 0);
      fCells.push_back(pv); fImp.push_back(std::pow(C.f, std::min(k + 1, n)));
    }
    fWorld = gw;
  }
  void CreateStore() {
    G4IStore* st = G4IStore::GetInstance(GetName());
    st->AddImportanceGeometryCell(1.0, *fWorld, 0);
    for (size_t i = 0; i < fCells.size(); ++i) st->AddImportanceGeometryCell(fImp[i], *fCells[i], 0);
  }
  G4VPhysicalVolume* fWorld = nullptr;
  std::vector<G4VPhysicalVolume*> fCells;
  std::vector<double> fImp;
};

// ---------------------------------------------------------------- fuente
class Gen : public G4VUserPrimaryGeneratorAction {
 public:
  Gen() { fGun = new G4ParticleGun(1); fGun->SetParticleDefinition(G4Neutron::Definition()); }
  void GeneratePrimaries(G4Event* ev) override {
    const Src& s = C.src[std::min((size_t)(G4UniformRand() * C.src.size()), C.src.size() - 1)];
    double r = G4UniformRand(); size_t i = 1; while (i < C.cdf.size() - 1 && C.cdf[i] < r) ++i;
    double f = (r - C.cdf[i - 1]) / std::max(C.cdf[i] - C.cdf[i - 1], 1e-30);
    double E = std::exp(std::log(C.E[i - 1]) + f * std::log(C.E[i] / C.E[i - 1]));
    double ct = 2 * G4UniformRand() - 1, ph = CLHEP::twopi * G4UniformRand(), st = std::sqrt(1 - ct * ct);
    fGun->SetParticleEnergy(E * MeV);
    fGun->SetParticlePosition(s.pos);
    fGun->SetParticleMomentumDirection(G4ThreeVector(st * std::cos(ph), st * std::sin(ph), ct));
    fGun->GeneratePrimaryVertex(ev);
  }
 private:
  G4ParticleGun* fGun;
};

// ---------------------------------------------------------------- puntuación (por evento → media y varianza)
static std::vector<double> evN, evG, sN, sN2, sG, sG2;
class Step : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* st) override {
    const G4VPhysicalVolume* pv = st->GetPreStepPoint()->GetPhysicalVolume();
    if (pv == gVault && st->GetPreStepPoint()->GetPosition().y() > C.iy / 2) {   // diagnóstico por capa (muro frontal)
      int k = (int)((st->GetPreStepPoint()->GetPosition().y() - C.iy / 2) / 100.0);
      if (k >= 0 && k < 40) { gLayN[k] += 1; gLayW[k] += st->GetStepLength() / cm * st->GetPreStepPoint()->GetWeight(); }
    }
    if (pv == gVault) {                 // ruleta "dura": fuera del semiespacio del muro frontal (y < +iy/2), las
      const G4ThreeVector& x = st->GetPostStepPoint()->GetPosition();   // partículas a > 40 cm dentro del hormigón
      const double d = 400.0;                                            // no vuelven (albedo despreciable)
      if (x.y() < C.iy / 2 && (std::abs(x.x()) > C.ix / 2 + d || x.y() < -C.iy / 2 - d ||
                               x.z() > C.floor_z + C.ih + d || x.z() < C.floor_z - d))
        st->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }
    auto it = DIDX.find(pv);
    if (it == DIDX.end()) return;
    const G4Track* tr = st->GetTrack();
    double E = st->GetPreStepPoint()->GetKineticEnergy() / MeV;
    double L = st->GetStepLength() / cm * st->GetPreStepPoint()->GetWeight() / D[it->second].vol_cm3;
    if (tr->GetDefinition() == G4Neutron::Definition()) evN[it->second] += L * loglog(nE, nH, 47, E);
    else if (tr->GetDefinition() == G4Gamma::Definition()) evG[it->second] += L * loglog(gE, gH, 25, E);
  }
};
class Ev : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override { std::fill(evN.begin(), evN.end(), 0); std::fill(evG.begin(), evG.end(), 0); }
  void EndOfEventAction(const G4Event*) override {
    for (size_t i = 0; i < D.size(); ++i) { sN[i] += evN[i]; sN2[i] += evN[i] * evN[i]; sG[i] += evG[i]; sG2[i] += evG[i] * evG[i]; }
  }
};
// solo se transportan neutrones y gammas (los cargados depositan localmente: no cuentan para H*(10) externo)
class Stack : public G4UserStackingAction {
 public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
    auto* p = t->GetDefinition();
    return (p == G4Neutron::Definition() || p == G4Gamma::Definition()) ? fUrgent : fKill;
  }
};
class Act : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new Gen); SetUserAction(new Step); SetUserAction(new Ev); SetUserAction(new Stack);
  }
};

int main(int argc, char** argv) {
  if (argc < 7) { std::cerr << "uso: srg_bench config gdml t_cm N semilla salida [f]\n"; return 1; }
  readCfg(argv[1]);
  C.t = std::atof(argv[3]) * 10.0;
  long N = std::atol(argv[4]);
  G4Random::setTheSeed(std::atol(argv[5]));
  if (argc > 7) C.f = std::atof(argv[7]);
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* geo = new Geo(argv[2]);
  auto* imp = new Imp;
  geo->RegisterParallelWorld(imp);
  rm->SetUserInitialization(geo);
  G4PhysListFactory fac;
  G4VModularPhysicsList* pl = fac.GetReferencePhysList("QGSP_BIC_HP");
  pl->SetVerboseLevel(0);
  auto* gsN = new G4GeometrySampler("imp", "neutron"); gsN->SetParallel(true);
  auto* gsG = new G4GeometrySampler("imp", "gamma");   gsG->SetParallel(true);
  pl->RegisterPhysics(new G4ImportanceBiasing(gsN, "imp"));
  pl->RegisterPhysics(new G4ImportanceBiasing(gsG, "imp"));
  pl->RegisterPhysics(new G4ParallelWorldPhysics("imp"));
  rm->SetUserInitialization(pl);
  rm->SetUserInitialization(new Act);
  rm->Initialize();
  imp->CreateStore();
  size_t n = 6; evN.assign(n, 0); evG.assign(n, 0); sN.assign(n, 0); sN2.assign(n, 0); sG.assign(n, 0); sG2.assign(n, 0);
  rm->BeamOn((G4int)N);
  std::ofstream o(argv[6]);
  o << "# srg_bench  t_cm " << C.t / 10 << "  N " << N << "  f " << C.f << "  seed " << argv[5] << "\n";
  o << "# detector  Hn_media_pSv_por_n  suma_cuad_N  Hg_media_pSv_por_n  suma_cuad_G  x y z\n";
  for (size_t i = 0; i < D.size(); ++i)
    o << D[i].name << " " << std::setprecision(10) << sN[i] / N << " " << sN2[i] << " " << sG[i] / N << " " << sG2[i]
      << " " << D[i].c.x() << " " << D[i].c.y() << " " << D[i].c.z() << "\n";
  o << "# capas del muro frontal (10 cm): pasos sin peso | longitud de traza ponderada [cm] por evento\n";
  for (int k = 0; k < 40; ++k) if (gLayN[k] > 0) o << "# capa " << k << " " << gLayN[k] << " " << gLayW[k] / N << "\n";
  delete rm;
  return 0;
}
