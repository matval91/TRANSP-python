import pytransp.classes.transp_exp as te
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels
common_style()

colors = ['r','b','g', 'c', 'm']
shot=58823;
data={
      '1':{'fname':'/home/vallar/TCV/58823/58823V69.CDF',\
           'time':[0.7, 1.]},
      '2':{'fname':'/home/vallar/TCV/58823/58823V74.CDF',\
           'time':[1, 1.15]},\
      '3':{'fname':'/home/vallar/TCV/58823/58823V73.CDF',\
           'time':[1.15, 1.35]},\
      '4':{'fname':'/home/vallar/TCV/58823/58823V71.CDF',\
           'time':[1.35, 2.]},\
}
fnames = [data[i]['fname'] for i in data.keys()]
labels = [r'$\tau=10$ ms, $z_{eff}$=1.8',\
          r'$\tau=19$ ms, $z_{eff}$=2',\
              #r'D=0.5 $m^2/s$, $z_{eff}$=3',\
          ]
exp_fname_w = '/home/vallar/TCV/58823/58823wmhd.dat'
exp_fname_v = '/home/vallar/TCV/58823/58823vl.dat'
exp_fname_n0 = '/home/vallar/TCV/58823/58823n0.dat'

wmhddata = np.loadtxt(exp_fname_w)
n0data = np.loadtxt(exp_fname_n0)

#Filtering on the two times where the wmhd matches

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(wmhddata[:,0], wmhddata[:,1]*1e-3, 'k', label='Exp.', lw=2.3)

f1=plt.figure()
axn0=f1.add_subplot(111)
axn0.plot(n0data[:,0], n0data[:,1], 'k', label='Exp.', lw=2.3)

wtot=np.array([])
n0=np.array([])
time=np.array([])
for i, key in enumerate(['1','2','3','4']):
    o=te.transp_exp(data[key]['fname'])
    if key=='1':
        tth, wth = o.t, o.wth
    o._calculate_wmhd()
    timet=data[key]['time']
    it = np.logical_and(o.t<timet[1], o.t>timet[0])
    time = np.append(time, o.t[it])                     
    wtot = np.append(wtot, o.wtot[it])
    o._calculate_n0()
    n0 = np.append(n0, o.n0_tot[it, -1])
ax.plot(time, wtot*1e-3, 'r', lw=2.3, label='Simulation, total')
ax.plot(tth, wth*1e-3, 'r--', lw=2.3, label='Simulation, thermal')
axn0.plot(time, n0*1e6, 'r', lw=2.3, label='Simulation')
    
limit_labels(ax, r'Time (s)', r'Energy [kJ]')
ax.legend(loc='best', fontsize='medium')
limit_labels(axn0, r'Time (s)',r'n$_0^{edge}$ [1/m$^3$]')
axn0.legend(loc='best', fontsize='medium')
f1.tight_layout()

ax.grid('on')
plt.show()
