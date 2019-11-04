#!/usr/bin/env python

import pytransp.plot.plot_eq as pe
import sys
import numpy as np
import pytransp.classes.transp_output as to
import matplotlib
import matplotlib.pyplot as plt
import import_tcv_eq as ite
matplotlib.interactive(True)

if len(sys.argv)==1:
    print("Please give as input the shots id and requested time. The simulations will be luckily found in /home/vallar/tr_client/TCV/<<shot>>")
    print('\n e.g. \n compare_transp_eq.py 62855V01 1.3 \n')
    sys.exit()

fnames = ['']*(len(sys.argv)-2)
_pre = "/home/vallar/NoTivoli/tr_client/TCV/"
pre = _pre+sys.argv[1][0:5]+"/"
suff=''
if sys.argv[1][-1]!="F":
    suff=".CDF"
fnames[0] = pre+sys.argv[1]+suff
time=float(sys.argv[-1])
shot=sys.argv[1][0:5]

print 'Reading files ', fnames
print 'At time ', time, 's' 

output=to.transp_output(fnames[0])

s=ite.get_tcv_eq(int(shot), time)
output.plot_equilibrium(time, 0)
f=plt.gcf()
rho=ite.psi_to_phi(s)
ite.plot_tcv_eq(s,f)
matplotlib.pyplot.show(block=True)
