#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 17:46:18 2019

@author: vallar
"""
import pytransp.classes.transp_fbm as tfbm
import sys
import matplotlib.pyplot as plt
col = ['r', 'g', 'b']
def main(fnames=['../examples/65052V01_fi_1.cdf']):
    """
    """
    for i,el in enumerate(fnames):
        fbm=tfbm.transp_fbm(el)
        fbm.plot_space()
        fbm.plot_Epitch()
        if i==0:
            fbm.plot_spaceE(); ax_pitch = plt.gca()
            fbm.plot_spacep(); ax_E = plt.gca()
        else:
            fbm.plot_spaceE(ax=ax_pitch, color=col[i]); 
            fbm.plot_spacep(ax=ax_E, color=col[i]);   
    return

if __name__ == '__main__':
    fnames=sys.argv[1:]
    main(fnames)