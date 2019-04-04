"""
Script to merge the two simulations with the different anomalous diffusion coefficient
"""
import pytransp.classes.transp_heating as th
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels, _cumulative_plot, _plot_1d
common_style()
col = ['k','r','b','m','g','c']
col2=np.copy(col)
col2[0]='y'
data={
      '1':{'fname':'/home/vallar/TCV/58823/58823V69.CDF',\
           'time':[0.7, 1.]},
      '2':{'fname':'/home/vallar/TCV/58823/58823V74.CDF',\
           'time':[1.01, 1.15]},\
      '3':{'fname':'/home/vallar/TCV/58823/58823V73.CDF',\
           'time':[1.16, 1.34]},\
      '4':{'fname':'/home/vallar/TCV/58823/58823V71.CDF',\
           'time':[1.35, 2.]},\
}
fnames = [data[i]['fname'] for i in data.keys()]

t=np.array([]); 
yabs=np.array([]); 
yfabs=np.array([]); 
ycx=np.array([]); 
yfcx=np.array([]); 
yol=np.array([]); 
yfol=np.array([]); 
yst=np.array([]); 
yfst=np.array([]); 
pe=np.array([]); pi=np.array([]); g=np.array([])
nbcd=np.array([]); unbcd=np.array([]); eff=np.array([])
peprof=np.zeros((2,50))
piprof=np.zeros((2,50))
nrel=np.zeros((2,50))
jprof=np.zeros((2,50))
rho=np.zeros((2,50))

for i, key in enumerate(['1','2','3','4']):
    o=th.transp_heating(data[key]['fname'])
    data[key]['tobj']=o
    tt=data[key]['time']
    it = np.logical_and(o.t<tt[1], o.t>tt[0])
    it=np.where(it==True)[0]
    t=np.append(t, o.t[it])
    x,y,yH,yfrac = o._init_pbalance_vars()
    yabs=np.append(yabs, y[0,it])
    ycx=np.append(ycx, y[1,it])
    yol=np.append(yol, y[2,it])
    yst=np.append(yst, y[3,it])
    yfabs=np.append(yfabs, yfrac[0,it])
    yfcx=np.append(yfcx, yfrac[1,it])
    yfol=np.append(yfol, yfrac[2,it])
    yfst=np.append(yfst, yfrac[3,it])
        
    pe=np.append(pe, o.pe[it]); pi=np.append(pi, o.pi[it]); 
    g=np.append(g, o.gi_vavg[it])
    nbcd = np.append(nbcd, o.nbcd[it])
    unbcd = np.append(unbcd, o.unbcd[it])
    eff = np.append(eff, o.eff[it])
   
    
    #choosing slices for plotting profiles
    for i,ttt in enumerate([0.9, 1.25]):
        if min(o.t[it])<ttt and max(o.t[it])>=ttt:
            indt = np.argmin(o.t-ttt  < 0.)
            rho[i,:] = o.rho[indt,:]
            peprof[i,:]=o.nb_FKP_vars['pe'][indt,:]
            piprof[i,:]=o.nb_FKP_vars['pi'][indt,:]
            nrel[i,:]=o.nb_FKP_vars['n'][indt,:]/o.kin_vars['ne'][indt,:]
            jprof[i,:]=o.nb_FKP_vars['nbcd'][indt,:]

ypbalance=np.array([yabs, ycx, yol, yst])
yfbalance=np.array([yfabs, yfcx, yfol, yfst]) 
###################
# Power balance
###################
labels = ['Abs.', 'CX', 'O. L.', 'Shine-thr.']
xlabel = r'Time [s]'
ylabel = r'Power [MW]'

f=plt.figure(figsize=(15,5))
axabs = f.add_subplot(121)
_cumulative_plot(t,ypbalance,labels, xlabel, ylabel, col2, ax=axabs, title='')

axfrac = f.add_subplot(122)
#axfrac = f.add_subplot(111)

for i in range(4):
    axfrac.plot(t, yfbalance[i,:]*100., color=col2[i], lw=2.3, label=labels[i])
axfrac.set_ylim([0,60])
axfrac.legend(loc='best', fontsize='medium')
#axabs.grid('on')
axfrac.grid('on')
f.subplots_adjust(right=0.8)        
f.tight_layout()

###################
# Power gi, ge
###################
f = plt.figure(figsize=(10, 8))
axp = f.add_subplot(111)
f=plt.figure()
axf = f.add_subplot(111, sharex=axp)
_plot_1d(t, pe*1e-3, ax=axp, color='k', label=r'e')
_plot_1d(t, pi*1e-3, ax=axp, color='r', label=r'i')
_plot_1d(t, pe/(pe+pi)*100., ax=axf, color='k', label=r'e')
_plot_1d(t, pi/(pe+pi)*100., ax=axf, color='r', label=r'i')
axf.plot(t, g*1e2, 'c', lw=2.3, label=r'G')

limit_labels(axp,r'Time [s]',r'P[kW]','')
limit_labels(axf,r'Time [s]',r'%','')
axp.set_ylim([0,220]); axf.set_ylim([30., 65.])
axf.legend(loc='best')
f.tight_layout()


###################
# Power deposition profiles (pe, pi, nf/ntot)
###################

f=plt.figure(figsize=(18,6)); 
axi = f.add_subplot(131); axe=f.add_subplot(132); axn=f.add_subplot(133)
labs=[r't=0.9 s',r't=1.25 s']; col=['k','r']
for i in [0,1]:
    lab=labs[i]
    axi.plot(rho[i,:], piprof[i,:]*1e-3,col[i], lw=2, label=lab)
    axe.plot(rho[i,:], peprof[i,:]*1e-3,col[i], lw=2, label=lab)
    axn.plot(rho[i,:], nrel[i, :]*100., col[i], lw=2, label=lab)

limit_labels(axe,r'$\rho$', r'$P_e$ [$kW/m^3$]','')
limit_labels(axi,r'$\rho$', r'$P_i$ [$kW/m^3$]','')
limit_labels(axn,r'$\rho$', r'$n_f/n_e$ [%]','')
#limit_labels(axj,r'$\rho$', r'j shielded [$kA/m^2$]','')
axn.legend(bbox_to_anchor=(1.05, 0.65), loc=2)
f.tight_layout()
f.subplots_adjust(right=0.8)
#
#
f = plt.figure(figsize=(10,8))
axj  = f.add_subplot(224)
axcd = f.add_subplot(221)
axsh = f.add_subplot(222)
axeff = f.add_subplot(223)
        
# for ax in [axcd, axsh, axeff]:
#     if time[0]!=0:
#         for i, el in enumerate(ind):
#             ax.axvline(x=self.t[self.inj_index[el]], color=col[i], lw=1.8, linestyle='--')

axcd.plot(t, nbcd*1e-3, 'k', lw=2.3)
axcd.plot(t, unbcd*1e-3, 'k--', lw=2.3)

axsh.plot(t, 1.-nbcd/unbcd, 'k', lw=2.3)
axeff.plot(t, eff*100., 'k', lw=2.3)
for i in [0,1]:
    lab=labs[i]
    axj.plot(rho[i,:], jprof[i,:]*1e-3,col[i], lw=2, label=lab)

axj.set_xlabel(r'$\rho$'); axj.set_ylabel('$j^{SH}$ [$kA/m^2$]')
limit_labels(axcd, r't [s]', r'$I_{CD}$ [$kA$]'); axj.grid('on')
axcd.legend(loc='best')
limit_labels(axsh, r't [s]', r'$1-I_{SH}/I_{UN}$')
limit_labels(axeff, r't [s]', r' $\eta \left[\frac{10^{18} A}{W m^2}\right]$', r'NBCD efficiency')
axeff.set_ylim([0,1.])
axj.legend(loc='upper right', fontsize='medium')
#f.text(0.01, 0.01, self.runid)
f.tight_layout();
#
plt.show()
