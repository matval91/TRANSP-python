import transp_output
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from a4py.utils.ascot_utils import common_style, limit_labels
common_style()

colors = ['r','b','g', 'c', 'm']

shot=58823;
#shot=58832
if shot==58823:
    fname1 = '/home/vallar/TCV/58823/58823V42.CDF' #TAUP=19ms D=0, my zeff
    fname2 = '/home/vallar/TCV/58823/58823V63.CDF' #TAUP=19ms D=0.5, my zeff
    # fname3 = '/home/vallar/TCV/58823/58823V52.CDF' #TAUP=19ms D=0.5, 2.0 constant
    # fname4 = '/home/vallar/TCV/58823/58823V45.CDF' 
    # fname5 = '/home/vallar/TCV/58823/58823V48.CDF'

    #fnames = [fname1, fname2, fname3, fname4]#, fname5]
    labels = [r'$\tau=10$ ms, $z_{eff}$=1.8',\
              r'$\tau=19$ ms, $z_{eff}$=2',\
              #r'D=0.5 $m^2/s$, $z_{eff}$=3',\
              ]
    fnames = [fname1, fname2]#, fname3, fname4, fname5]
    exp_fname = '/home/vallar/TCV/58823/58823n0.dat'

n0data = np.loadtxt(exp_fname)
f=plt.figure()
f1d = plt.figure(); ax1d=f1d.add_subplot(111)
#ft = plt.figure(); axt=ft.add_subplot(111)

ax=f.add_subplot(111)
ax.plot(n0data[:,0], n0data[:,1]*1e-16, 'k', label='EXP', lw=2.3)
t=0.99
time1=1.05
time1=0.8
for i, fname in enumerate(fnames):
    o=transp_output.output_1d(fname)
    o._n0edge()
    if '42' in fname:
        it = np.logical_and(o.t<time1, o.t>0.8)
        it = o.t>time1
        n0 = o.n0_tot[it, -1]
        time = o.t[it]
        ind = np.argmin(o.t-t+0.1 < 0.)
        label=labels[0]
    elif '63' in fname:
        it = np.logical_and(o.t>time1, o.t>0.8)
        n0 = o.n0_tot[it, -1]
        time = o.t[o.t>time1]
        ind = np.argmin(o.t-t-0.1 < 0.)
        label=labels[1]
    elif '43' in fname:
        it = np.logical_and(it, o.t>0.8)
        n0 = o.n0_tot[it, -1]
        time = o.t[it]
        ind = np.argmin(o.t-t-0.1 < 0.)
        label=labels[2]
    else:
        continue
                  
    ax.plot(time, n0*1e6*1e-16, colors[i], label=label, lw=2.3)
    ax1d.semilogy(o.rho, o.n0_tot[ind,:]*1e6*1e-16, colors[i], label=label+'| t={:.2f} s'.format(o.t[ind]), lw=2.3)
    ax1d.semilogy([1., 1.1], [2, 2], colors[i])
    #axt.plot(o.rho, o.file.variables['TAUPI'][50,:]*1e3, colors[i], label=fname[-12:-4])


limit_labels(ax,r'Time (s)',r'n$_0^{edge} [10^{16}/m^3]$')
ax.legend(loc='best', fontsize='medium')
ax.set_ylim([0, 5]); ax.set_xlim([0, 2.])
f.tight_layout(); ax.grid('on');

#limit_labels(ax1d,r'$\rho$',r'n$_0 [10^{16}/m^3$]')
ax1d.set_xlabel(r'$\rho$'); ax1d.set_ylabel(r'n$_0 [10^{16}/m^3$]')
ax1d.legend(loc='best', fontsize='medium')
ax1d.yaxis.set_major_formatter(ScalarFormatter())
f1d.tight_layout(); ax1d.grid('on')

# axt.set_xlabel(r'$\rho$'); axt.set_ylabel(r'$\tau$ [ms]')
# axt.legend(loc='best')
# ft.tight_layout(); axt.grid('on')

plt.show()
