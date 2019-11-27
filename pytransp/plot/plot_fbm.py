# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pytransp.classes.transp_fbm as tfbm
import sys

    
def main(fname='../examples/65052V01_fi_1.cdf'):
    """
    """
    fbm=tfbm.transp_fbm(fname)
    fbm.plot_space()
    fbm.plot_Epitch()
    fbm.plot_spaceE()
    fbm.plot_spacep()
    return fbm

if __name__ == '__main__':
    fname=sys.argv[1]
    main(fname)