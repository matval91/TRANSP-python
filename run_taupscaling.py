#!/usr/bin/env python

import pytransp.compare_wexp.scale_taup as st
import sys
import numpy as np
import glob
if len(sys.argv)==1:
    print("Please give as input the shot id. The simulations will be searched in /home/vallar/tr_client/TCV/<<shot>> /n The n0 data in /home/vallar/edge_density /n     to produce the n0 data you can use save_n0 function inside /home/vallar/matlabscripts")
    
    print('\n e.g. \n run_taupscaling.py 62124 \n')
    sys.exit()

shot=sys.argv[1]
folder='/home/vallar/tr_client/TCV/'+shot+'/'
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

exp_fname='/home/vallar/edge_density/'+shot+'n0.dat'
try:
    f=open(fname)
    f.close()
except NameError:
    print('Files do not exist in '+folder)
    sys.exit()

st.scale_taup(fname, exp_fname, time=1., setvalue=1e25)
