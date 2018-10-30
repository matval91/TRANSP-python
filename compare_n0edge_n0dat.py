import transp_output
import numpy as np
import matplotlib.pyplot as plt

fname1 = '/home/vallar/TCV/58823/58823O03.CDF'
fname2 = '/home/vallar/TCV/58823/58823O05.CDF'
fname3 = '/home/vallar/TCV/58823/58823O06.CDF' #TAUP=5ms
fnames = [fname1, fname2, fname3]
colors = ['r','b','g']
exp_fname = '/home/vallar/TCV/58823/58823n0.dat'

n0data = np.loadtxt(exp_fname)
f=plt.figure()
f1d = plt.figure(); ax1d=f1d.add_subplot(111)
ax=f.add_subplot(111)
ax.plot(n0data[:,0], n0data[:,1], 'k--', label='EXP')

for i, fname in enumerate(fnames):
    o=transp_output.output_1d(fname)
    o._n0edge()
    ax.plot(o.t, o.n0_tot[:,-1]*1e6, colors[i], label=fname[-12:-4])
    ax1d.plot(o.rho, o.n0_tot[50,:]*1e6, colors[i], label=fname[-12:-4])


ax.set_xlabel(r'Time (s)'); ax.set_ylabel(r'n$_0^{edge}$ [1/m$^3$]')
ax.legend(loc='best')
f.tight_layout()

ax1d.set_xlabel(r'$\rho$'); ax1d.set_ylabel(r'n$_0$ [1/m$^3$]')
ax1d.legend(loc='best')
f1d.tight_layout()
ax.grid('on'); ax1d.grid('on')
plt.show()
