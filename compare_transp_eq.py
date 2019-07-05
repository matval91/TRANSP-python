#!/usr/bin/env python

import pytransp.plot.plot_eq as pe
import sys
import numpy as np
import pytransp.classes.transp_output as to
import matplotlib

matplotlib.interactive(True)

if len(sys.argv)==1:
    print("Please give as input the shots id and requested time. The simulations will be luckily found in /home/vallar/tr_client/TCV/<<shot>>")
    print('\n e.g. \n compare_transp_eq.py 62855V01 62855V02 1.3 \n')
    sys.exit()

fnames = ['']*(len(sys.argv)-2)
for i, el in enumerate(sys.argv[1:-1]):
    _pre = "/home/vallar/tr_client/TCV/"
    pre = _pre+el[0:5]+"/"
    suff=''
    if el[-1]!="F":
        suff=".CDF"
    fnames[i] = pre+el+suff
time=float(sys.argv[-1])

print 'Reading files ', fnames
print 'At time ', time, 's' 

output=np.empty((len(fnames)), dtype='object')
for i,el in enumerate(fnames):
    output[i]=to.transp_output(el)

for i,el in enumerate(output):
    if i==0: f=0
    else:
        f=plt.gcf()

    el.plot_equilibrium(time, f)
    
matplotlib.pyplot.show(block=True)
