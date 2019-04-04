import pytransp.classes.transp_exp as te
import numpy as np
import matplotlib.pyplot as plt

colors = ['r','b','g', 'c', 'm']
shot=58823
#shot=58832
if shot==58823:
    fname1 = '/home/vallar/TCV/58823/58823V71.CDF' #TAUP=5ms
    fname2 = '/home/vallar/TCV/58823/58823V65.CDF' #TAUP=5ms
    fname3 = '/home/vallar/TCV/58823/58823V66.CDF' #TAUP=5ms
    fname4 = '/home/vallar/TCV/58823/58823V70.CDF' #TAUP=5ms
    fname5 = '/home/vallar/TCV/58823/58823V69.CDF' #TAUP=5ms
    #fnames = [fname1, fname2,fname3, 
    fnames=[fname1]
    exp_fname = '/home/vallar/TCV/58823/58823n0.dat'
elif shot==58832:
    fname1 = '/home/vallar/TCV/58832/58832V19.CDF' #tauph=20ms
    fname2 = '/home/vallar/TCV/58832/58832V20.CDF'
    fname3 = '/home/vallar/TCV/58832/58832V21.CDF' #TAUP=5ms
    fname4 = '/home/vallar/TCV/58832/58832V22.CDF' #TAUP=5ms
    fname5 = '/home/vallar/TCV/58823/58823V23.CDF' #TAUP=5ms
    
    fnames = [fname1, fname2,fname3, fname4,fname5]
    exp_fname = '/home/vallar/TCV/58832/58832n0.dat'

n0data = np.loadtxt(exp_fname)
f=plt.figure()
f1d = plt.figure(); ax1d=f1d.add_subplot(111)
#ft = plt.figure(); axt=ft.add_subplot(111)

ax=f.add_subplot(111)
ax.plot(n0data[:,0], n0data[:,1], 'k--', label='EXP')
t=0.9
for i, fname in enumerate(fnames):
    try:
        o=te.transp_exp(fname)
    except:
        continue
    ind = np.argmin(o.t-t < 0.)
    o._calculate_n0()
    ax.plot(o.t, o.n0_tot[:,-1]*1e6, colors[i], label=fname[-12:-4])
    ax1d.semilogy(o.rho[ind,:], o.n0_tot[ind,:]*1e6, colors[i], label=fname[-12:-4])
    ax1d.semilogy([1., 1.1], [2e16, 2e16], colors[i])
    #axt.plot(o.rho, o.file.variables['TAUPI'][50,:]*1e3, colors[i], label=fname[-12:-4])


ax.set_xlabel(r'Time (s)'); ax.set_ylabel(r'n$_0^{edge}$ [1/m$^3$]')
ax.legend(loc='best')
ax.set_ylim([0, 5e16]); ax.set_xlim([0, 2.])
f.tight_layout(); ax.grid('on');

ax1d.set_xlabel(r'$\rho$'); ax1d.set_ylabel(r'n$_0$ [1/m$^3$]')
ax1d.legend(loc='best')
f1d.tight_layout(); ax1d.grid('on')

# axt.set_xlabel(r'$\rho$'); axt.set_ylabel(r'$\tau$ [ms]')
# axt.legend(loc='best')
# ft.tight_layout(); axt.grid('on')

plt.show()
