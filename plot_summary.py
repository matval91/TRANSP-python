#!/usr/bin/env python


import pytransp.plot.summary_plot as sp
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.integrate as integrate
import glob, sys

if len(sys.argv) >= 2:
    shot = sys.argv[1]
    if len(sys.argv) == 3:
        tid = sys.argv[2]
else:
    print("Please give as input shot number, id (OPTIONAL)")
    print("Files will be looked for in '/home/vallar/tr_client/TCV/<<SHOT>>")
    print('\n e.g. \n plot_summary.py 62124 (V01) \n')
    sys.exit()

folder='/home/vallar/NoTivoli/tr_client/TCV/'+str(shot)+'/'
try:
    tid
    fname=folder+str(shot)+str(tid)+'.CDF'
except NameError:
    files = [f for f in glob.glob(folder + "*.CDF")]
    _np=0
    for i,f in enumerate(files):
        _n = f[-6:-4]
        if int(_n)<_np:
            continue
        else:
            _np=int(_n)
            ind_f=i
    try:
        ind_f
    except NameError:
        print('Files do not exist in '+folder)
        sys.exit()
        
    fname=files[ind_f]        
print('Getting CDF file '+fname)

sp.summary_plot(fname)
