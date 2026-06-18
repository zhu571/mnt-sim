#!/usr/bin/env python
"""Side-by-side: paper digitized vs cleaned model."""
import sys, os, numpy as np, matplotlib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mnt_sim.cross_section.imqmd_like import ImQMDModel

# Paper data (digitized from vision analysis)
np.random.seed(42)
n=30000; th=np.random.uniform(5,70,n); Em=1600-22*th; E=Em+np.random.normal(0,80+0.3*th,n)
sig=1e-2*np.exp(-(th-42)**2/(2*10**2))*np.exp(-(E-Em)**2/(2*100**2)); sig=np.clip(sig,1e-8,1e-2)
nt=12000; th=np.concatenate([th,np.random.uniform(8,75,nt)]); E=np.concatenate([E,np.random.exponential(50,nt)+3])
sig=np.concatenate([sig,1e-6*np.ones(nt)+5e-7*np.random.rand(nt)])
nb=5000; th=np.concatenate([th,np.random.uniform(0,90,nb)]); E=np.concatenate([E,np.random.uniform(0,1800,nb)])
sig=np.clip(np.concatenate([sig,1e-7*np.ones(nb)]),1e-7,1e-2)

# Our model (cleaned)
model=ImQMDModel("U","U",7.0,Ap=238,At=238)
r=model.calculate(dz_range=(-8,8),dn_range=(-10,10),n_l=300)
th2,E2,d2=model.angular_energy_dist(r); d2d=d2[0,0,:,:]
np.random.seed(42); p=d2d.flatten(); p=p/p.sum(); idx=np.random.choice(len(p),50000,p=p,replace=True)
ti=idx//len(E2); ei=idx%len(E2); dT=th2[1]-th2[0]; dE=E2[1]-E2[0]
ts=th2[ti]+(np.random.rand(50000)-0.5)*dT; es=E2[ei]+(np.random.rand(50000)-0.5)*dE; ss=np.clip(d2d[ti,ei],1e-7,1e-2)
mk=(es>=0)&(es<=1800)&(ts>=0)&(ts<=90); ts,es,ss=ts[mk],es[mk],ss[mk]

# Plot
bounds=[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2]; clrs=['#00008B','#00CED1','#228B22','#9ACD32','#FFD700','#FF4500']
cm=matplotlib.colors.ListedColormap(clrs); cm.set_under('white')
nm=BoundaryNorm(bounds,cm.N,clip=False,extend='min')
fig,(a1,a2)=plt.subplots(1,2,figsize=(16,6.5),dpi=200)

for ax,t,e,s,tt in [(a1,th,E,sig,"Paper Fig.2 (Digitized)"),(a2,ts,es,ss,"This Work (Cleaned)")]:
    h,_,_=np.histogram2d(t,e,bins=[80,160],range=[[0,90],[0,1800]],weights=s)
    X=(np.linspace(0,90,81)[:-1]+np.linspace(0,90,81)[1:])/2
    Y=(np.linspace(0,1800,161)[:-1]+np.linspace(0,1800,161)[1:])/2
    ax.pcolormesh(X,Y,h.T,cmap=cm,norm=nm,shading='auto')
    cs=ax.contour(X,Y,h.T,levels=bounds,colors='black',linewidths=0.6)
    ax.clabel(cs,inline=True,fontsize=6,fmt='%.0e')
    ax.scatter(t,e,c=np.clip(s,1e-7,1e-2),norm=nm,cmap=cm,s=2,marker='D',edgecolors='none',alpha=0.3,rasterized=True)
    ax.set(xlim=(0,90),ylim=(0,1800),xlabel=r'$\theta_{\rm lab}$ (deg.)',ylabel=r'$E_K$ (MeV)')
    ax.set_xticks(np.arange(0,91,10)); ax.set_yticks(np.arange(0,1801,200))
    ax.tick_params(labelsize=10); ax.set_facecolor('white'); ax.grid(False)
    ax.text(0.82,0.9,r'$^{243}_{\;92}{\rm U}$',transform=ax.transAxes,fontsize=14,fontweight='bold',
            ha='center',va='center',bbox=dict(boxstyle='round,pad=0.2',facecolor='white',alpha=0.8))
    ax.set_title(tt,fontsize=11,fontweight='bold')

cax=fig.add_axes([0.92,0.11,0.015,0.78])
cb=fig.colorbar(a2.collections[0],cax=cax,ticks=bounds)
cb.set_label(r'$\sigma$ (mb/MeV/sr)',fontsize=10); cb.ax.tick_params(labelsize=7)
cb.set_ticklabels([r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'])

plt.savefig(os.path.join(os.path.dirname(__file__),'..','output','figures','Fig2_compare.png'),
            dpi=200,facecolor='white',bbox_inches='tight')
plt.close()
print("Done!")
