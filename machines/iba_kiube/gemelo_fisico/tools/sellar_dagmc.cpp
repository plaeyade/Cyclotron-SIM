// sellar_dagmc6.cpp — topología DAGMC completa del búnker BIND:
// volumen pared (mat:Concreto) + graveyard (mat:Graveyard) + cavidad como
// complemento implícito + boundary:vacuum en la caja exterior.
// Convención GEOM_SENSE_2 = {forward, reverse}; las facetas apuntan hacia
// fuera del volumen forward.
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include <iostream>
#include <cstring>
#include <cmath>
using namespace moab;

static void chk(ErrorCode r, const char* paso) {
  if (r != MB_SUCCESS) { std::cerr << "FALLO " << paso << " rc=" << r << "\n"; std::exit(9); }
}

int main(int argc, char** argv) {
  if (argc < 3) { std::cerr << "uso: sellar6 in.h5m out.h5m\n"; return 1; }
  Core mb;
  chk(mb.load_file(argv[1]), "load");

  Tag geom, idt, sense, name, cat;
  chk(mb.tag_get_handle("GEOM_DIMENSION", 1, MB_TYPE_INTEGER, geom), "geom");
  chk(mb.tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, idt), "id");
  chk(mb.tag_get_handle("GEOM_SENSE_2", 2, MB_TYPE_HANDLE, sense,
                        MB_TAG_SPARSE | MB_TAG_CREAT), "sense");
  chk(mb.tag_get_handle("NAME", 32, MB_TYPE_OPAQUE, name,
                        MB_TAG_BYTES|MB_TAG_SPARSE|MB_TAG_CREAT, nullptr, nullptr), "name");
  chk(mb.tag_get_handle("CATEGORY", 32, MB_TYPE_OPAQUE, cat,
                        MB_TAG_BYTES|MB_TAG_SPARSE|MB_TAG_CREAT, nullptr, nullptr), "cat");

  // bounding box del modelo (muros) en las unidades del fichero (m)
  double bbmin[3]={1e30,1e30,1e30}, bbmax[3]={-1e30,-1e30,-1e30};
  { Range vs; mb.get_entities_by_type(0, MBVERTEX, vs);
    for (auto v : vs) { double c[3]; mb.get_coords(&v,1,c);
      for (int a=0;a<3;a++){ if(c[a]<bbmin[a])bbmin[a]=c[a]; if(c[a]>bbmax[a])bbmax[a]=c[a]; } } }
  double cen[3]; for (int a=0;a<3;a++) cen[a]=(bbmin[a]+bbmax[a])/2;
  double umbral = 0.30*std::min(bbmax[0]-bbmin[0], bbmax[1]-bbmin[1]);

  Range surfs;
  const int two=2; const void* v2[]={&two};
  chk(mb.get_entities_by_type_and_tag(0, MBENTITYSET, &geom, v2, 1, surfs), "surfs");
  std::cout << "superficies muro: " << surfs.size() << "\n";

  auto tagset = [&](EntityHandle h, int dim, int gid, const char* nm, const char* ct) {
    chk(mb.tag_set_data(geom,&h,1,&dim), "t_dim");
    chk(mb.tag_set_data(idt,&h,1,&gid), "t_id");
    char b1[32]={0}; std::strncpy(b1,nm,31); chk(mb.tag_set_data(name,&h,1,b1), "t_nm");
    char b2[32]={0}; std::strncpy(b2,ct,31); chk(mb.tag_set_data(cat,&h,1,b2), "t_ct");
  };

  // volumen pared (id 1) y graveyard (id 2)
  EntityHandle wall, grave;
  chk(mb.create_meshset(MESHSET_SET, wall), "wall");
  chk(mb.create_meshset(MESHSET_SET, grave), "grave");
  tagset(wall, 3, 1, "wall", "Volume");
  tagset(grave, 3, 2, "graveyard", "Volume");

  // caja graveyard: 1 m más allá del muro en cada eje (unidades del fichero)
  double gx0=bbmin[0]-1.0, gx1=bbmax[0]+1.0;
  double gy0=bbmin[1]-1.0, gy1=bbmax[1]+1.0;
  double gz0=bbmin[2]-1.0, gz1=bbmax[2]+1.0;
  double vc[8][3]={{gx0,gy0,gz0},{gx1,gy0,gz0},{gx1,gy1,gz0},{gx0,gy1,gz0},
                   {gx0,gy0,gz1},{gx1,gy0,gz1},{gx1,gy1,gz1},{gx0,gy1,gz1}};
  EntityHandle vv[8];
  for (int i=0;i<8;i++) chk(mb.create_vertex(vc[i], vv[i]), "vtx");
  // 6 caras (normales hacia fuera, regla mano derecha), 2 triángulos c/u
  int faces[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
  int next_surf_id = 100;
  Range gsurfs;
  for (int f=0; f<6; f++) {
    EntityHandle s; chk(mb.create_meshset(MESHSET_SET, s), "gsurf");
    tagset(s, 2, next_surf_id++, "gsurf", "Surface");
    EntityHandle t1, t2;
    const EntityHandle q1[3]={vv[faces[f][0]],vv[faces[f][1]],vv[faces[f][2]]};
    const EntityHandle q2[3]={vv[faces[f][0]],vv[faces[f][2]],vv[faces[f][3]]};
    chk(mb.create_element(MBTRI, q1, 3, t1), "tri1");
    chk(mb.create_element(MBTRI, q2, 3, t2), "tri2");
    chk(mb.add_entities(s, &t1, 1), "a1"); chk(mb.add_entities(s, &t2, 1), "a2");
    // sentido: facetas hacia fuera de la caja -> graveyard en reverse
    EntityHandle sv[2]={0, grave};
    chk(mb.tag_set_data(sense, &s, 1, sv), "gsense");
    chk(mb.add_parent_child(grave, s), "gpc");
    gsurfs.insert(s);
  }

  // clasificar caras del muro: exteriores (compartidas con graveyard) vs cavidad
  Range outer, cavity;
  for (auto s : surfs) {
    Range tris; mb.get_entities_by_type(s, MBTRI, tris);
    double cx=0,cy=0,cz=0; int n=0;
    for (auto t : tris) { const EntityHandle* c; int nn; mb.get_connectivity(t,c,nn);
      for (int k=0;k<nn;k++){ double p[3]; mb.get_coords(&c[k],1,p); cx+=p[0];cy+=p[1];cz+=p[2];n++; } }
    cx/=n; cy/=n; cz/=n;
    // cara exterior si su centroide esta pegado al bbox en su eje dominante
    double ex=cx-cen[0], ey=cy-cen[1], ez=cz-cen[2];
    double ax=std::fabs(ex), ay=std::fabs(ey), az=std::fabs(ez);
    double tol=0.05;
    bool ext;
    if (ax>=ay && ax>=az) ext = (cx < bbmin[0]+tol || cx > bbmax[0]-tol);
    else if (ay>=ax && ay>=az) ext = (cy < bbmin[1]+tol || cy > bbmax[1]-tol);
    else ext = (cz < bbmin[2]+tol || cz > bbmax[2]-tol);
    if (ext) outer.insert(s); else cavity.insert(s);
  }
  std::cout << "exteriores: " << outer.size() << " cavidad: " << cavity.size() << "\n";

  int sid=1;
  for (auto s : cavity) {  // normales hacia el muro: {0, wall}
    EntityHandle sv[2]={0, wall};
    chk(mb.tag_set_data(sense,&s,1,&sv), "s_cav");
    chk(mb.tag_set_data(idt,&s,1,&sid), "id_cav");
    char ct[32]={0}; std::strncpy(ct,"Surface",31); chk(mb.tag_set_data(cat,&s,1,ct),"c_cav");
    chk(mb.add_parent_child(wall, s), "pc_cav"); sid++;
  }
  for (auto s : outer) {   // normales hacia dentro (hacia cavidad): grave forward, wall reverse
    EntityHandle sv[2]={grave, wall};
    chk(mb.tag_set_data(sense,&s,1,&sv), "s_out");
    chk(mb.tag_set_data(idt,&s,1,&sid), "id_out");
    char ct[32]={0}; std::strncpy(ct,"Surface",31); chk(mb.tag_set_data(cat,&s,1,ct),"c_out");
    chk(mb.add_parent_child(wall, s), "pc_w");
    chk(mb.add_parent_child(grave, s), "pc_g"); sid++;
  }

  // grupos UWUW
  auto grupo = [&](int gid, const char* nm, EntityHandle miembro, bool is_vol) {
    EntityHandle g; chk(mb.create_meshset(MESHSET_SET, g), "grp");
    tagset(g, 4, gid, nm, "Group");
    chk(mb.add_entities(g, &miembro, 1), "grp_add");
    return g;
  };
  grupo(1, "mat:Concreto/rho:2.3", wall, true);
  grupo(2, "mat:Graveyard", grave, true);
  // boundary:vacuum sobre las caras exteriores de la caja graveyard
  EntityHandle bc; chk(mb.create_meshset(MESHSET_SET, bc), "bc");
  tagset(bc, 4, 3, "boundary:vacuum", "Group");
  for (auto s : gsurfs) chk(mb.add_entities(bc, &s, 1), "bc_add");

  // tolerancia de facetado
  Tag ftol; chk(mb.tag_get_handle("FACETING_TOL", 1, MB_TYPE_DOUBLE, ftol,
                                  MB_TAG_SPARSE|MB_TAG_CREAT), "ftol");
  const double ft = 0.001;
  Range sets; mb.get_entities_by_type(0, MBENTITYSET, sets);
  for (auto st : sets) mb.tag_set_data(ftol, &st, 1, &ft);

  chk(mb.write_file(argv[2]), "write");
  std::cout << "ESCRITO " << argv[2] << "\n";
  return 0;
}
