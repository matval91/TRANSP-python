import transp_output
import numpy as np
import matplotlib.pyplot as plt

fname1 = '/home/vallar/TCV/58823/58823V45.CDF' #TAUP=5ms
fname2 = '/home/vallar/TCV/58823/58823V41.CDF' #TAUP=5ms
fname1 = '/home/vallar/TCV/58823/58823V37.CDF' #TAUP=5ms
fname2 = '/home/vallar/TCV/58823/58823V38.CDF' #TAUP=5ms
fname3 = '/home/vallar/TCV/58823/58823V42.CDF' #TAUP=5ms
fname4 = '/home/vallar/TCV/58823/58823V43.CDF' #TAUP=5ms
fname5 = '/home/vallar/TCV/58823/58823V44.CDF' #TAUP=5ms

fnames = [fname1, fname2, fname3, fname4, fname5]
#fnames = [fname1, fname5]
colors = ['r','b','g', 'c', 'm']

f=plt.figure(); ax=f.add_subplot(111)
fn=plt.figure(); axn=fn.add_subplot(111)
fi=plt.figure(); axi=fi.add_subplot(111)
time=1.25
for i, fname in enumerate(fnames):
    o=transp_output.output_1d(fname)
    ind = np.argmin(o.t - time < 0.)
    print(fname); print('zeff',o.file.variables['ZEFFC'][ind])
    rho = o.rho; pbe = o.nb_FKP_vars['pe'][ind,:]; n = o.nb_FKP_vars['n'][ind,:]
    pbi = o.nb_FKP_vars['pi'][ind,:];    
    ax.plot(rho, pbe, colors[i], label=fname[-12:-4])
    axi.plot(rho, pbi, colors[i], label=fname[-12:-4])
    axn.plot(rho, n*1e6, colors[i], label=fname[-12:-4])

ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'pbe')
ax.legend(loc='best')
f.tight_layout(); ax.grid('on');

axn.set_xlabel(r'$\rho$'); axn.set_ylabel(r'n [1/m$^3$]')
axn.legend(loc='best')
fn.tight_layout(); axn.grid('on');

axi.set_xlabel(r'$\rho$'); axi.set_ylabel(r'pbi')
axi.legend(loc='best')
fi.tight_layout(); axi.grid('on');

plt.show()
