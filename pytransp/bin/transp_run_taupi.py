#!/usr/bin/env python

import pytransp.compare_wexp.scale_taup as st
import sys, os
import numpy as np
import glob
import ufiles as uf
import matplotlib.pyplot as plt

flag_run=False
if len(sys.argv)==1:
    print("Please give as input the shot id. The simulations will be searched in /home/vallar/tr_client/TCV/<<shot>> /n The n0 data in /home/vallar/edge_density /n     to produce the n0 data you can use save_n0 function inside /home/vallar/matlabscripts")
    
    print('\n e.g. \n run_taupscaling.py 62124 \n')
    sys.exit()
elif len(sys.argv)>2:
    flag_run=True

shot=sys.argv[1]
folder='/home/vallar/NoTivoli/tr_client/TCV/'+shot+'/'
files = [f for f in glob.glob(folder + "*.CDF")]

if flag_run==False:
    _np=0
    for i,f in enumerate(files):
        _n = f[-6:-4]
        if int(_n)<_np:
            continue
        else:
            _np=int(_n)
            ind_f=i
else:
    ind_f=1

try:
    ind_f
except NameError:
    print('Files do not exist in '+folder)
    sys.exit()

fname=files[ind_f]        
print('Getting CDF file '+fname)

## Produce edge density file from baratrons
os.system('matlab -nodisplay -r "cd(\'{}\'); save_n0({});exit"'.format('/home/vallar/edge_density', shot))
command = 'cp /tmp/vallar/'+shot+'n0.dat /home/vallar/edge_density/'+shot+'n0.dat'
print('Doing command '+command)
os.system(command)
exp_fname='/home/vallar/edge_density/'+shot+'n0.dat'
try:
    f=open(fname)
    f.close()
except NameError:
    print('Files do not exist in '+folder)
    sys.exit()

st.scale_taup(fname, exp_fname, time=1., setvalue=1e25)

tpi=uf.RU('/home/vallar/OMF58499.TPI')
f=plt.figure(); ax=f.add_subplot(111)
ax.plot(tpi.values['X'], tpi.fvalues*1e-3, lw=2.3)
ax.set_xlabel(r't (s)'); ax.set_ylabel(r'$\tau_i$ [ms]')
