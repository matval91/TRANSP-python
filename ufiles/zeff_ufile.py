"""
Script to produce Zeff for TRANSP starting from a mat file
"""

import scipy.io as sio
import ufiles as uf
import sys

try:
    shot = sys.argv[1]
except:
    shot = '58782'

print "Shot number: "+str(shot)
fname = '/home/vallar/MST1/T16/icdb/zeffin_'+shot+'.mat'
struct_name = fname[-16:-4]
zeff_s = sio.loadmat(fname)
data = zeff_s[struct_name][0,0]
time = data[0][:,0]
zeff = data[1][:,0]

uf_d = {'pre': 'OMF','ext': 'ZFF',
        'shot': shot,
        'grid': {'X': {'lbl': ' TIME                SEC       ',
                       'arr': time}},
        'data': {'lbl':  'Effective Charge    None ',
                 'arr': zeff}}
 
uf.WU(uf_d, udir='/home/vallar')
