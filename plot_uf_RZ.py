"""
Script to plot R and Z surf from ufiles R*****.SURF and Z*****.SURF
"""
from __future__ import print_function
import ufiles as uf
import numpy as np
import matplotlib.pyplot as plt

dir = '/home/vallar/tr_client/TCV/'
shot = 59331
shot = 58823
dir += str(shot)+'/'

Rs=uf.RU(dir+'R'+str(shot)+'.SURF')
zs=uf.RU(dir+'Z'+str(shot)+'.SURF')
R   = Rs.fvalues*0.01
z   = zs.fvalues*0.01
tf  = Rs.values['X']
the = Rs.values['Y']
rho = Rs.values['Z']

time=1.01
tind = np.argmin(np.abs(time-tf)>0)
print(tf[tind])
f=plt.figure(figsize=(5, 8))
ax = f.add_subplot(111)
ax.set_title(str(shot)+' t='+str(time))
ax.set_xlabel('R [m]')
ax.set_ylabel('Z [m]')
for i,rr in enumerate(rho):
    ax.plot(R[tind, :, i], z[tind, :, i], lw=1.5)

ax.axis('equal')
f.tight_layout()
plt.show()
