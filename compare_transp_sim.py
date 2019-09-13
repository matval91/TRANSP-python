#!/usr/bin/env python

import pytransp.plot.summary_plot as sp
import sys
import numpy as np
import matplotlib
from utils.plot_utils import common_style

common_style()
matplotlib.interactive(True)
matplotlib.pyplot.rc('xtick', labelsize=15)
matplotlib.pyplot.rc('ytick', labelsize=15)
matplotlib.pyplot.rc('legend', fontsize=15)

if len(sys.argv)==1:
    print("Please give as input the shots id and requested time. The simulations will be searched in /home/vallar/tr_client/TCV/<<shot>>")
    print('\n e.g. \n compare_transp_sim.py 62855V01 62855V02 1.3 \n')
    sys.exit()

fnames = ['']*(len(sys.argv)-2)
for i, el in enumerate(sys.argv[1:-1]):
    _pre = "/home/vallar/NoTivoli/tr_client/TCV/"
    pre = _pre+el[0:5]+"/"
    suff=''
    if el[-1]!="F":
        suff=".CDF"
    fnames[i] = pre+el+suff
time=float(sys.argv[-1])

print 'Reading files ', fnames
print 'At time ', time, 's' 

sp.summary_plot(fname_sim=fnames, t=time)
matplotlib.pyplot.show(block=True)
