#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Script to compare two fbm distributions

Created on Wed Jan  8 15:44:22 2020

@author: vallar
"""
import pytransp.classes.transp_fbm as tf
import matplotlib.pyplot as plt
import utils.plot_utils as pu
import sys
import numpy as np

def main(fname1, fname2, R, z):
    f1=tf.transp_fbm(fname1)
    f2=tf.transp_fbm(fname2)
    f1.plot_Epitch(); f2.plot_Epitch()
    #mean=f1.f_space_int
    x=f1.dict_dim['pitch']; y=f2.dict_dim['E']*1e-3
    indR = np.argmin(f1.dict_dim['R']-R<0.)
    indz = np.argmin(f1.dict_dim['z']-z<0.)    
    
    diff=f2.fdist_notnorm[:,:,indR,indz]-f1.fdist_notnorm[:,:,indR,indz]
    
    f=plt.figure()
    ax=f.add_subplot(111)
    cb=ax.contourf(x,y,diff,  cmap='rainbow')
    ax.set_xlabel(r'$\xi$')
    ax.set_ylabel(r'E [keV]')
    plt.colorbar(cb)
    pu.limit_labels(ax, r'$\xi$',r'E [keV]', title='(t2-t1,R,z)=({:.2f}-{:.2f} s,{:.2f} m,{:.2f} m)'.format(f2.time, f1.time, R,z))
    f.tight_layout()
    ax.grid('on')
    plt.show()

if len(sys.argv) ==5:
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    R = sys.argv[3]
    z = sys.argv[4]
    main(fname1, fname2, float(R), float(z))
else:
    print("Please give as input the time of interest and a/some CDF shot files")
    print("")
    print('\n e.g. \n plot_nubeam.py 1. 65052V01.CDF \n')
    sys.exit()
    
    
