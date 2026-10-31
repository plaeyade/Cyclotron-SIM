// bunker_iba (RPOS F89): búnker del ciclotrón con la GEOMETRÍA IBA del PDF GAMMA (bóveda 4.0×3.94×2.5 m, muros O 2.40 / E 2.20 /
// S 2.50, techo 2.30, losa 0.55, laberinto esquemático de 5 tramos al norte), concreto de alta densidad estándar 2.5 g/cm³,
// 300 µA en UN blanco (T1, orientado hacia la radiofarmacia = este). Importancia PLANAR en una cara por corrida (este, oeste,
// sur, techo) → perfiles a través de cada barrera. Coordenadas mundo: este=+y, oeste=-y, sur=+x, norte=-x, techo=+z.
// Derivado de bunker_mapa v1b — equivalente v0.35 de RPOS F70: máquina KIUBE (GDML del LOD neutrónico, v0.35) en el búnker del proyecto
// (bloque BIND de hormigón: muros, losa de piso y techo como un solo casco), fuente 18O(p,n) desde los blancos F18
// T1 y T5 en la cara del yugo con las puertas de HDPE (BON21, isótropa), física QGSP_BIC_HP.
// Reducción de varianza: importancia geométrica en cascarones laterales anidados de 10 cm (factor f por cascarón) a
// través de los muros → hay cuentas fuera de los muros (F70 tenía 0 y se extrapoló en F78). Techo y losa sin sesgo.
// Tallies: (1) malla 3D de H*(10) n y γ (ICRP 74, longitud de traza con recorrido exacto de vóxeles); (2) fluencia de
// neutrones integrada por región (hormigón por capa de profundidad y por cara; acero; cobre de bobinas; cobre de dees;
// inox) en grupos térmico (ponderado 1/v), epitérmico (0.5 eV–0.1 MeV) y rápido (> 1 MeV) para la activación.
// Uso: bunker_mapa config.txt kiube.gdml N semilla salida_prefijo [f=2]
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

#include "G4Region.hh"
// ---------------------------------------------------------------- configuración
struct Src { std::string name; G4ThreeVector pos, dir; };
struct Cfg {
  double floor_z = -1120, ix = 5900, iy = 4700, ih = 2500, wall = 2110, slab = 550, roof = 1450, rho = 2.35, f = 2.0;
  double m0[3] = {-6060, -5460, -2670}, m1[3] = {6060, 5460, 3830}, mh = 100;   // malla [mm]
  std::vector<std::pair<std::string, double>> wt;
  std::vector<Src> src;
  std::vector<double> E, cdf;
};
static Cfg C;
static void readCfg(const char* fn) {
  std::ifstream in(fn); std::string line; std::vector<double> w;
  while (std::getline(in, line)) {
    std::istringstream s(line); std::string k; s >> k;
    if (k == "floor_z") s >> C.floor_z;
    else if (k == "vault_inner") s >> C.ix >> C.iy >> C.ih;
    else if (k == "wall") s >> C.wall;
    else if (k == "slab") s >> C.slab;
    else if (k == "roof") s >> C.roof;
    else if (k == "mesh") s >> C.m0[0] >> C.m1[0] >> C.m0[1] >> C.m1[1] >> C.m0[2] >> C.m1[2] >> C.mh;
    else if (k == "concrete_density") s >> C.rho;
    else if (k == "concrete_wt") { std::string e; double v; while (s >> e >> v) C.wt.push_back({e, v}); }
    else if (k == "source") { Src x; double a, b, c, d, e; s >> x.name >> a >> b >> c >> d >> e;
      x.pos = G4ThreeVector(a, b, c); x.dir = G4ThreeVector(d, e, 0); C.src.push_back(x); }
    else if (k == "spectrum") { int n; s >> n;
      for (int i = 0; i < n; ++i) { std::getline(in, line); std::istringstream t(line); double e, y; t >> e >> y; C.E.push_back(e); w.push_back(y); } }
  }
  C.cdf.assign(C.E.size(), 0.0);
  for (size_t i = 1; i < C.E.size(); ++i) C.cdf[i] = C.cdf[i - 1] + 0.5 * (w[i] + w[i - 1]) * std::log(C.E[i] / C.E[i - 1]);
  for (auto& c : C.cdf) c /= C.cdf.back();
}
// distancia (mm) de un punto a la sala interior (0 dentro); cara: 0 -x(oeste) 1 +x(este) 2 -y(sur) 3 +y(norte) 4 piso 5 techo
static double depth(const G4ThreeVector& x, int& face) {
  double d[6] = {-x.x() - C.ix / 2, x.x() - C.ix / 2, -x.y() - C.iy / 2, x.y() - C.iy / 2, C.floor_z - x.z(), x.z() - (C.floor_z + C.ih)};
  face = 0; for (int i = 1; i < 6; ++i) if (d[i] > d[face]) face = i;
  return std::max(0.0, d[face]);
}

// ---------------------------------------------------------------- malla y tallies
static int NX, NY, NZ; static size_t NV;
static std::vector<float> mN, mG;            // suma por corrida (pSv·cm por n fuente, sin dividir por volumen)
static std::vector<float> fT, fE, fF;         // v1b: fluencia de neutrones (cm por n fuente): <0.5 eV, 0.5 eV–0.1 MeV, >0.1 MeV
static G4VPhysicalVolume* gBunker = nullptr;
static std::map<const G4Material*, int> MREG; static std::vector<std::string> RNAME;
static const int NLAY = 30;                    // capas de 10 cm del hormigón
// por región: térmico 1/v, epitérmico, rápido, total [cm por n fuente]; hormigón por (cara, capa)
static std::vector<double> rTh, rEp, rFa, rTo;
static double cTh[6][NLAY], cEp[6][NLAY], cFa[6][NLAY], cTo[6][NLAY];

static void meshAdd(const G4ThreeVector& a, const G4ThreeVector& b, double w, std::vector<float>& M) {
  // recorrido exacto del segmento a→b por la malla (Amanatides–Woo); w = valor por cm
  G4ThreeVector d = b - a; double L = d.mag(); if (L <= 0) return;
  double t0 = 0, t1 = 1;
  for (int k = 0; k < 3; ++k) {                // recorte a la caja de la malla
    double p = a[k], q = d[k];
    if (std::abs(q) < 1e-12) { if (p < C.m0[k] || p > C.m1[k]) return; continue; }
    double ta = (C.m0[k] - p) / q, tb = (C.m1[k] - p) / q; if (ta > tb) std::swap(ta, tb);
    t0 = std::max(t0, ta); t1 = std::min(t1, tb); if (t0 >= t1) return;
  }
  G4ThreeVector p = a + t0 * d; int n[3] = {NX, NY, NZ}, i[3], st[3]; double tm[3], td[3];
  for (int k = 0; k < 3; ++k) {
    i[k] = std::min(n[k] - 1, std::max(0, (int)std::floor((p[k] - C.m0[k]) / C.mh)));
    if (d[k] > 0) { st[k] = 1; tm[k] = (C.m0[k] + (i[k] + 1) * C.mh - a[k]) / d[k]; td[k] = C.mh / d[k]; }
    else if (d[k] < 0) { st[k] = -1; tm[k] = (C.m0[k] + i[k] * C.mh - a[k]) / d[k]; td[k] = -C.mh / d[k]; }
    else { st[k] = 0; tm[k] = 1e30; td[k] = 1e30; }
  }
  double t = t0;
  while (t < t1) {
    int k = (tm[0] < tm[1]) ? (tm[0] < tm[2] ? 0 : 2) : (tm[1] < tm[2] ? 1 : 2);
    double tn = std::min(tm[k], t1);
    M[((size_t)i[0] * NY + i[1]) * NZ + i[2]] += (float)(w * (tn - t) * L / 10.0);
    t = tn; if (t >= t1) break;
    i[k] += st[k]; tm[k] += td[k];
    if (i[k] < 0 || i[k] >= n[k]) break;
  }
}

class Geo : public G4VUserDetectorConstruction {
 public:
  explicit Geo(const char* gdml) : fG(gdml) {}
  G4VPhysicalVolume* Construct() override {
    G4GDMLParser p; p.SetOverlapCheck(false); p.Read(fG, false);
    G4VPhysicalVolume* W = p.GetWorldVolume(); G4LogicalVolume* WL = W->GetLogicalVolume();
    auto* nist = G4NistManager::Instance();
    auto* conc = new G4Material("concreto_bind", C.rho * g / cm3, (G4int)C.wt.size());
    double tot = 0; for (auto& e : C.wt) tot += e.second;
    for (auto& e : C.wt) conc->AddElement(nist->FindOrBuildElement(e.first), e.second / tot);
    // planta IBA (m, origen en el centro de la bóveda): X_plan (este+) → y mundo; Y_plan (norte+) → −x mundo
    const double fz = C.floor_z, zb = fz - C.slab, zt = fz + C.ih + C.roof;
    auto wbox = [&](double X0, double X1, double Y0, double Y1, double z0, double z1, const char* n) {
      auto* b = new G4Box(n, (Y1 - Y0) * 500.0, (X1 - X0) * 500.0, (z1 - z0) / 2);
      G4ThreeVector c(-(Y0 + Y1) * 500.0, (X0 + X1) * 500.0, (z0 + z1) / 2); return std::make_pair(b, c); };
    auto out = wbox(-4.40, 4.20, -4.47, 9.21, zb, zt, "bo");
    G4VSolid* sh = out.first; G4ThreeVector co = out.second;
    struct R { double X0, X1, Y0, Y1; };
    std::vector<R> huecos = {{-2.0, 2.0, -1.97, 1.97},                       // bóveda
                             {0.9, 2.0, 1.96, 2.94},                          // entrada al laberinto
                             {-3.2, 2.0, 2.93, 3.93}, {-3.2, -2.2, 3.92, 4.44},   // tramo 1 + enlace oeste
                             {-3.2, 2.2, 4.43, 5.43}, {1.2, 2.2, 5.42, 5.94},     // tramo 2 + enlace este
                             {-3.2, 2.2, 5.93, 6.93}, {-3.2, -2.2, 6.92, 7.44},   // tramo 3 + enlace oeste
                             {-3.2, 3.2, 7.43, 8.43}, {2.2, 3.2, 8.42, 9.22}};    // tramo 4 + salida (puerta sin blindaje)
    int i = 0;
    for (auto& h : huecos) {
      auto b = wbox(h.X0, h.X1, h.Y0, h.Y1, fz, fz + C.ih, ("h" + std::to_string(i++)).c_str());
      sh = new G4SubtractionSolid("bunker" + std::to_string(i), sh, b.first, nullptr, b.second - co);
    }
    gBunker = new G4PVPlacement(nullptr, co, new G4LogicalVolume(sh, conc, "bunker"), "bunker", WL, false, 0, true);
    // regiones de activación por material de la máquina
    const char* mats[] = {"steel_1010", "copper_coil_homog", "copper_ofe", "copper_dee_caja", "stainless_304L", "hdpe", "niobium", "water_o18"};
    for (auto* m : *G4Material::GetMaterialTable()) {
      std::string nm = m->GetName();
      for (const char* k : mats) if (nm.find(k) != std::string::npos && !MREG.count(m)) {
        int idx = -1; for (size_t j = 0; j < RNAME.size(); ++j) if (RNAME[j] == k) idx = (int)j;
        if (idx < 0) { RNAME.push_back(k); idx = (int)RNAME.size() - 1; }
        MREG[m] = idx; }
    }
    return W;
  }
 private:
  const char* fG;
};

static int AX = 1; static double SG = 1, P0 = 2000, NP = 32;      // eje (0 x,1 y,2 z), signo, plano inicial [mm], nº planos
class Imp : public G4VUserParallelWorld {
  // importancia planar: celda k = franja entre los planos P0+100k y P0+100(k+1) a lo largo de SG·eje (resto del mundo
  // perpendicular completo); importancia f^(k+1); antes del plano 0 importancia 1
 public:
  Imp() : G4VUserParallelWorld("imp") {}
  void Construct() override {
    G4VPhysicalVolume* gw = GetWorld(); G4LogicalVolume* gl = gw->GetLogicalVolume();
    auto* wb = (G4Box*)gl->GetSolid();
    double H[3] = {wb->GetXHalfLength() - 1, wb->GetYHalfLength() - 1, wb->GetZHalfLength() - 1};
    for (int k = 0; k < (int)NP; ++k) {
      double h[3] = {H[0], H[1], H[2]}; h[AX] = 50.0;
      G4ThreeVector c(0, 0, 0); c[AX] = SG * (P0 + 100.0 * k + 50.0);
      auto* L = new G4LogicalVolume(new G4Box("p" + std::to_string(k), h[0], h[1], h[2]), nullptr, "p" + std::to_string(k));
      fCells.push_back(new G4PVPlacement(nullptr, c, L, "p" + std::to_string(k), gl, false, 0));
    }
    fWorld = gw;
  }
  void CreateStore() {
    G4IStore* st = G4IStore::GetInstance(GetName());
    st->AddImportanceGeometryCell(1.0, *fWorld, 0);
    for (size_t i = 0; i < fCells.size(); ++i) st->AddImportanceGeometryCell(std::pow(C.f, (double)(i + 1)), *fCells[i], 0);
  }
  G4VPhysicalVolume* fWorld = nullptr; std::vector<G4VPhysicalVolume*> fCells;
};

class Gen : public G4VUserPrimaryGeneratorAction {
 public:
  Gen() { fGun = new G4ParticleGun(1); fGun->SetParticleDefinition(G4Neutron::Definition()); }
  void GeneratePrimaries(G4Event* ev) override {
    const Src& s = C.src[std::min(C.src.size() - 1, (size_t)(G4UniformRand() * C.src.size()))];   // F90: doble haz = blancos equiprobables
    double r = G4UniformRand(); size_t i = 1; while (i < C.cdf.size() - 1 && C.cdf[i] < r) ++i;
    double f = (r - C.cdf[i - 1]) / std::max(C.cdf[i] - C.cdf[i - 1], 1e-30);
    double E = std::exp(std::log(C.E[i - 1]) + f * std::log(C.E[i] / C.E[i - 1]));
    double ct = 2 * G4UniformRand() - 1, ph = CLHEP::twopi * G4UniformRand(), st = std::sqrt(1 - ct * ct);
    fGun->SetParticleEnergy(E * MeV); fGun->SetParticlePosition(s.pos);
    fGun->SetParticleMomentumDirection(G4ThreeVector(st * std::cos(ph), st * std::sin(ph), ct));
    fGun->GeneratePrimaryVertex(ev);
  }
 private:
  G4ParticleGun* fGun;
};

class Step : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* st) override {
    const G4Track* tr = st->GetTrack();
    const bool isN = tr->GetDefinition() == G4Neutron::Definition(), isG = tr->GetDefinition() == G4Gamma::Definition();
    if (!isN && !isG) return;
    const G4StepPoint* a = st->GetPreStepPoint();
    double E = a->GetKineticEnergy() / MeV, w = a->GetWeight(), Lcm = st->GetStepLength() / cm;
    if (Lcm <= 0) return;
    double h = isN ? loglog(nE, nH, 47, E) : loglog(gE, gH, 25, E);
    meshAdd(a->GetPosition(), st->GetPostStepPoint()->GetPosition(), w * h, isN ? mN : mG);
    if (!isN) return;
    meshAdd(a->GetPosition(), st->GetPostStepPoint()->GetPosition(), w, E < 5e-7 ? fT : (E < 0.1 ? fE : fF));
    double th = (E < 5e-7) ? Lcm * w * std::sqrt(2.53e-8 / std::max(E, 1e-12)) : 0.0;   // 1/v (σ0-equivalente)
    double ep = (E >= 5e-7 && E < 0.1) ? Lcm * w : 0.0, fa = (E >= 1.0) ? Lcm * w : 0.0, to = Lcm * w;
    const G4VPhysicalVolume* pv = a->GetPhysicalVolume();
    if (pv == gBunker) {
      int face; double d = depth(a->GetPosition(), face); int k = std::min(NLAY - 1, (int)(d / 100.0));
      cTh[face][k] += th; cEp[face][k] += ep; cFa[face][k] += fa; cTo[face][k] += to; return;
    }
    auto it = MREG.find(a->GetMaterial());
    if (it != MREG.end()) { int r = it->second; rTh[r] += th; rEp[r] += ep; rFa[r] += fa; rTo[r] += to; }
  }
};
class Stack : public G4UserStackingAction {
 public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
    auto* p = t->GetDefinition();
    return (p == G4Neutron::Definition() || p == G4Gamma::Definition()) ? fUrgent : fKill;
  }
};
class Act : public G4VUserActionInitialization {
 public:
  void Build() const override { SetUserAction(new Gen); SetUserAction(new Step); SetUserAction(new Stack); }
};

int main(int argc, char** argv) {
  if (argc < 6) { std::cerr << "uso: bunker_iba config gdml N semilla salida_prefijo f eje signo plano0_mm nplanos\n"; return 1; }
  readCfg(argv[1]);
  long N = std::atol(argv[3]); G4Random::setTheSeed(std::atol(argv[4]));
  if (argc > 6) C.f = std::atof(argv[6]);
  if (argc > 9) { AX = std::atoi(argv[7]); SG = std::atof(argv[8]); P0 = std::atof(argv[9]); NP = (argc > 10) ? std::atof(argv[10]) : 32; }
  NX = (int)std::lround((C.m1[0] - C.m0[0]) / C.mh); NY = (int)std::lround((C.m1[1] - C.m0[1]) / C.mh); NZ = (int)std::lround((C.m1[2] - C.m0[2]) / C.mh);
  NV = (size_t)NX * NY * NZ; mN.assign(NV, 0.f); mG.assign(NV, 0.f); fT.assign(NV, 0.f); fE.assign(NV, 0.f); fF.assign(NV, 0.f);
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* geo = new Geo(argv[2]); auto* imp = new Imp; geo->RegisterParallelWorld(imp);
  rm->SetUserInitialization(geo);
  G4PhysListFactory fac; G4VModularPhysicsList* pl = fac.GetReferencePhysList("QGSP_BIC_HP"); pl->SetVerboseLevel(0);
  auto* gsN = new G4GeometrySampler("imp", "neutron"); gsN->SetParallel(true);
  auto* gsG = new G4GeometrySampler("imp", "gamma"); gsG->SetParallel(true);
  pl->RegisterPhysics(new G4ImportanceBiasing(gsN, "imp")); pl->RegisterPhysics(new G4ImportanceBiasing(gsG, "imp"));
  pl->RegisterPhysics(new G4ParallelWorldPhysics("imp"));
  rm->SetUserInitialization(pl); rm->SetUserInitialization(new Act);
  rm->Initialize(); imp->CreateStore();
  rTh.assign(RNAME.size(), 0); rEp.assign(RNAME.size(), 0); rFa.assign(RNAME.size(), 0); rTo.assign(RNAME.size(), 0);
  rm->BeamOn((G4int)N);
  std::string pre = argv[5];
  { std::ofstream o(pre + "_malla.bin", std::ios::binary);
    int hdr[3] = {NX, NY, NZ}; o.write((char*)hdr, sizeof hdr);
    double g[7] = {C.m0[0], C.m0[1], C.m0[2], C.mh, (double)N, 0, 0}; o.write((char*)g, sizeof g);
    o.write((char*)mN.data(), NV * sizeof(float)); o.write((char*)mG.data(), NV * sizeof(float));
    o.write((char*)fT.data(), NV * sizeof(float)); o.write((char*)fE.data(), NV * sizeof(float)); o.write((char*)fF.data(), NV * sizeof(float)); }
  std::ofstream o(pre + "_regiones.txt");
  o << "# bunker_mapa N " << N << " f " << C.f << " seed " << argv[4] << "  [cm de traza por n fuente: térmico 1/v, epitérmico 0.5eV-0.1MeV, rápido >1MeV, total]\n";
  for (size_t r = 0; r < RNAME.size(); ++r)
    o << "maquina " << RNAME[r] << " " << std::setprecision(10) << rTh[r] / N << " " << rEp[r] / N << " " << rFa[r] / N << " " << rTo[r] / N << "\n";
  const char* fn[6] = {"oeste", "este", "sur", "norte", "piso", "techo"};
  for (int f = 0; f < 6; ++f) for (int k = 0; k < NLAY; ++k)
    if (cTo[f][k] > 0) o << "hormigon " << fn[f] << " " << k << " " << std::setprecision(10) << cTh[f][k] / N << " " << cEp[f][k] / N << " " << cFa[f][k] / N << " " << cTo[f][k] / N << "\n";
  delete rm;
  return 0;
}
