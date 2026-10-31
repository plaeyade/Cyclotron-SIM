"""Gap efectivo de la cámara de vacío (v0.35): volumen de vacío entre caras de polo (|z| ≤ 100 mm, r < 548 mm) / π r².
Comparación con el supuesto KIUBE-V04 (gap 132 mm) del repositorio RPOS. Salida: validation/gap_efectivo_camara.json"""
import sys, math, json; sys.path.insert(0,'.')
from kiube import load_params
from kiube.beam import extraction
from kiube.mc import mc_volumes
from kiube.geometry import zcyl
P,_=load_params(); E=extraction(P)
vols,_=mc_volumes(P,E,"openmc")
vac=next(v for v in vols if v["mat"]=="vacuum")["shape"]
R=548.0
cyl=zcyl(R,-100.001,100.001)
inner=vac.intersect(cyl)
V=inner.Volume(); bb=inner.BoundingBox()
A=math.pi*R*R
res=dict(V_vacio_total_m3=vac.Volume()/1e9, V_camara_r548_m3=V/1e9, gap_efectivo_mm=V/A, z_min=bb.zmin, z_max=bb.zmax)
# por anillos
for r0,r1 in ((60,200),(200,400),(400,548)):
    ring=inner.intersect(zcyl(r1,-100.001,100.001).cut(zcyl(r0,-100.001,100.001)))
    res[f"gap_ef_{r0}_{r1}_mm"]=ring.Volume()/(math.pi*(r1*r1-r0*r0))
print(json.dumps(res))
json.dump(res,open('validation/gap_efectivo_camara.json','w'), indent=1)
