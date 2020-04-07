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


def main(fname1, fname2):
    f1=tf.transp_fbm(fname1)
    f2=tf.transp_fbm(fname2)
    f1.plot_Epitch(); f2.plot_Epitch()
    diff=f1.f_space_int-f2.f_space_int
    #mean=f1.f_space_int
    x=f1.dict_dim['pitch']; y=f2.dict_dim['E']*1e-3
    
    
    f=plt.figure()
    ax=f.add_subplot(111)
    cb=ax.contourf(x,y,diff,  cmap='rainbow')
    ax.set_xlabel(r'$\xi$')
    ax.set_ylabel(r'E [keV]')
    plt.colorbar(cb)
    pu.limit_labels(ax, r'$\xi$',r'E [keV]', title='Post crash - pre crash')
    f.tight_layout()
    ax.grid('on')
    plt.show()

if len(sys.argv) ==3:
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    main(fname1, fname2)
else:
    print("Please give as input the time of interest and a/some CDF shot files")
    print("")
    print('\n e.g. \n plot_nubeam.py 1. 65052V01.CDF \n')
    sys.exit()
    
    
