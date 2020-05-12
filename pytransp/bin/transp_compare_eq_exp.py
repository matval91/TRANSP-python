#!/usr/bin/env python

import pytransp.plot.plot_eq as pe
import sys, glob
import numpy as np
import pytransp.classes.transp_output as to
import matplotlib
import matplotlib.pyplot as plt
import import_tcv_eq as ite
matplotlib.interactive(True)
liuqe_flavour= 'LIUQE.M'

if len(sys.argv) >= 3:
    shot = sys.argv[1]
    time = float(sys.argv[2])
    if len(sys.argv) == 4:
        tid = sys.argv[3]
    if len(sys.argv)== 5:
        liuqe_flavour= sys.argv[4]
else:
    print("Please give as input the shot and requested time. The simulations will be luckily found in /home/vallar/tr_client/TCV/<<shot>>\n")
    print('LIUQE_FLAVOUR = LIUQE.M, LIUQE.M2, LIUQE.M3 \n')
    print('\n e.g. \n compare_transp_eq.py 62855 1.3  (ID) (LIUQE_FLAVOUR)\n')
    sys.exit()

folder='/home/vallar//NoTivoli/tr_client/TCV/'+str(shot)+'/'
try:
    tid
    fname=folder+str(shot)+str(tid)+'.CDF'
except NameError:
    files = [f for f in glob.glob(folder + "*.CDF")]
    _np=0; ind_f=0;
    for i,f in enumerate(files):
        _n = f[-6:-4]
        if int(_n)<_np:
            continue
        else:
            _np=int(_n)
        ind_f=i
    fname=files[ind_f]

print('Getting CDF file '+fname)
print 'At time ', time, 's' 

output=to.transp_output(fname)

s=ite.get_tcv_eq(int(shot), time, liuqe_flavour=liuqe_flavour)
output.plot_equilibrium(time, 0)
rho=ite.psi_to_phi(s)

f=plt.gcf()
ite.plot_tcv_eq(s,f)
matplotlib.pyplot.show(block=True)
