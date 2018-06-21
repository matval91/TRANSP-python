"""
Script to interpolate Ti bad data from OMF ufile
input: shot number, initial time, end time
there must be a folder ~/omfit_<<SHOT>>/ with file /OMF<<SHOT>>.TI2
The data will be interpolated (linearly) between tini, tend.
"""
import scipy.interpolate as interp
import ufiles as uf
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    shot = sys.argv[1]
except:
    shot = '58782'
t = np.array([float(sys.argv[2]), float(sys.argv[3])])
    
print "Shot number: "+str(shot)
fname = '/home/vallar/omfit_'+shot+'/OMF'+shot+'.TI2'
uu = uf.OMFRU(fname)
rho = uu.values['X0']
t_uf = uu.values['X1']


ind_t = np.where(np.logical_and(t_uf>=t[0], t_uf<=t[1]))[0]
xt = np.array([t_uf[ind_t[0]], t_uf[ind_t[-1]]])
yt = np.array([uu.fvalues[:,ind_t[0]], uu.fvalues[:,ind_t[-1]]])
yt = np.transpose(yt)
param = interp.interp2d(xt, uu.values['X0'], yt)
xint = t_uf[ind_t]
yint = param(xint, rho)
f=plt.figure()
axold = f.add_subplot(211)
for iind, i in enumerate(ind_t):
    ls='--'
    if iind==0 or iind==len(ind_t)-1:
        ls='-'
    axold.plot(rho, uu.fvalues[:, i], linestyle=ls,label='{:.2f} s'.format(t_uf[i]))
plt.legend()
axnew = f.add_subplot(212)
for i in range(len(xint)):
    axnew.plot(rho, yint[:,i])
plt.show()

uu.fvalues[:, ind_t] = yint

uf_d = {'pre': 'OMF','ext': 'TI2',
         'shot': shot,
         'grid': {'X0': {'lbl': ' X      Poloidal Rho',
                         'arr': uu.values['X0']},
                  'X1': {'lbl': ' TIME               SEC       ',
                         'arr': uu.values['X1']}},
         'data': {'lbl':  'Ion Temperature   KEV ',
                  'arr': uu.fvalues}}

uf.OMFWU(uf_d, udir='/home/vallar')
