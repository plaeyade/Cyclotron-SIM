import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
R=json.load(open("out/mc_iba/resultado_iba.json")); C=R["caras"]; K=R["corridas"]
nom={"este":"Este (radiofarmacia)","oeste":"Oeste (exterior)","sur":"Sur, tapón (área verde)","techo":"Techo (área verde vacía)"}
fig,axs=plt.subplots(1,4,figsize=(16,4.6),sharey=True)
for ax,k in zip(axs,["este","oeste","sur","techo"]):
    r=C[k]; P=[p for p in r["perfil"] if p["prof_cm"]>0 and p["H"]>0]
    z=np.array([p["prof_cm"] for p in P]); H=np.array([p["H"] for p in P]); M=np.array([p["parche_max"] for p in P])
    ax.semilogy(z,H,"o-",ms=3,color="#1f4fd6",label="promedio en la cara")
    ax.semilogy(z,M,"-",lw=1,color="#d62728",alpha=.8,label="máx. parche 50×50 cm")
    ax.axhline(r["objetivo_uSv_h"],color="k",ls=":",lw=1); ax.text(3,r["objetivo_uSv_h"]*1.6,f"objetivo {r['objetivo_uSv_h']} µSv/h",fontsize=8)
    ax.axvline(r["espesor_IBA_cm"],color="#7f7f7f",lw=1.2); ax.text(r["espesor_IBA_cm"]+2,1e6,f"IBA {r['espesor_IBA_cm']:.0f}",fontsize=8,color="#555")
    m=r.get("modulo_13cm_directo")
    if m: ax.axvline(m,color="#2ca02c",lw=1.6); ax.text(m-4,1e4,f"req. {m} cm",fontsize=8,color="#2ca02c",ha="right")
    ax.axvspan(r["espesor_IBA_cm"],z.max()+10,color="#eeeeee",zorder=0)
    ax.set_title(f"{nom[k]}\n{K[k]['lotes']} lotes · {int(K[k]['historias'])} n fuente",fontsize=9.5)
    ax.set_xlabel("profundidad desde la cara interior [cm]"); ax.set_xlim(0,z.max()+10); ax.grid(alpha=.25,which="both",lw=.4)
axs[0].set_ylabel("H*(10) con haz [µSv/h], 300 µA"); axs[0].legend(fontsize=8,loc="lower left")
axs[0].set_ylim(1e-7,1e8)
fig.suptitle("Búnker ciclotrón IBA, tabique 2.5 (concreto estándar): Geant4 QGSP_BIC_HP, ¹⁸O(p,n), un blanco",fontsize=11,x=0.01,ha="left")
plt.tight_layout(); plt.savefig("out/mc_iba/perfiles_iba.png",dpi=160)
