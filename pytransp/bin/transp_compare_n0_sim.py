#!/usr/bin/env python

import pytransp.classes.transp_exp as te
import sys
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style

common_style()

if len(sys.argv)==1:
    print("Please give as input the shots id and requested time. The simulations will be searched in /home/vallar/tr_client/TCV/<<shot>>")
    print('\n e.g. \n compare_transp_n0.py 62855V01 62855V02 1.3 \n')
    sys.exit()

fnames = ['']*(len(sys.argv)-2)
   
for i, el in enumerate(sys.argv[1:-1]):
    _pre = "/home/vallar/tr_client/TCV/"
    try:
        output = te.transp_exp(_pre+el[0:5]+'/'+el+'.CDF');
    except IOError:
        _pre = "/home/vallar/NoTivoli/tr_client/TCV/"
    pre = _pre+el[0:5]+"/"
    suff=''
    if el[-1]!="F":
        suff=".CDF"
    fnames[i] = pre+el+suff
time=float(sys.argv[-1])

print 'Reading files ', fnames
print 'At time ', time, 's'

f=plt.figure();ax=f.add_subplot(111)
ax.set_xlabel(r'$\rho_{Tor}$'); ax.set_ylabel(r'$n_0$ [1/m$^3$]')
for el in fnames:
    output = te.transp_exp(el);
    ind = np.argmin(output.t-time<0.)
    x, y = output.rho[ind,:], output.n0_tot[ind,:]
    try:
        output._calculate_scalars()
        pcx = output.pcx[ind]/output.nb_in_vars['P_D'][ind]*100.
    except AttributeError:
        pcx=0.
    plt.plot(x, y, label=el[-12:-4]+' pcx='+'{:.2f}'.format(pcx)+' %')
plt.yscale('log')
plt.grid('on')
plt.legend(loc='best')
plt.tight_layout()

f=plt.figure();ax=f.add_subplot(111)
ax.set_xlabel(r'$\rho_{Tor}$'); ax.set_ylabel(r'$T_0^{rec}$ [eV]')
for el in fnames:
    output = te.transp_exp(el);
    ind = np.argmin(output.t-time<0.)
    x, y = output.rho[ind,:], output.file.variables['T0CX_RCD'][ind,:]
    try:
        output._calculate_scalars()
        pcx = output.pcx[ind]/output.nb_in_vars['P_D'][ind]*100.
    except AttributeError:
        pcx=0.
    plt.plot(x, y, label=el[-12:-4]+' pcx='+'{:.2f}'.format(pcx)+' %')
plt.yscale('linear')
plt.grid('on')
plt.legend(loc='best')
plt.tight_layout()


plt.show()
