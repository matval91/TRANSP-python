#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 13:31:23 2018

@author: vallar
"""

import pytransp.classes.transp_exp as te
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import ufiles as uf

colors = ['r','b','g', 'c', 'm']

def scale_taup(fname='/home/vallar/TCV/58823/58823V68.CDF', \
               exp_fname='/home/vallar/TCV/58823/58823n0.dat', \
               time=1., setvalue=1e25):
    """scales taup
    
    Produces a OMF***.TPI ufile with the confinement time as function of time.
    It tries to match the n0edge in exp_fname
    
    Parameters:
        fname (str): name of the TRANSP simulation
        exp_fname (str): name of the edge density (time, n0edge (m^-3))
        time (float) : if needed, time after where set the n0edge to a fixed value
        setvalue(float): value to set the edge neutral density
    Returns:
        None
    """

    taup=input('taup for simulation '+fname+' in millisec ')
    taup=float(taup)
    
    n0data = np.loadtxt(exp_fname)
    _texp, _n0exp = n0data[:,0], n0data[:,1]
    if setvalue!=1e25:
        _n0exp[_texp>time] = np.full(len(np.where(_texp>time)[0]), setvalue) 

    o=te.transp_exp(fname)
    o._calculate_n0()
    tsim, n0sim = o.t, o.n0_tot[:,-1]*1e6

    plt.plot(tsim, n0sim)
    plt.plot(_texp, _n0exp)
    

    param_exp = interp.interp1d(_texp, _n0exp)
    time = tsim[tsim<2.]
    n0exp = param_exp(time)
    

    factor = n0exp/n0sim
    taunew = taup/factor*1e-3 #converts to seconds

    f=plt.figure(); ax=f.add_subplot(111)
    ax.plot(time, n0exp, 'k--', label='EXP')
    ax.plot(time, n0sim, 'r-', label='SIM')
    ax.legend(loc='best')
    f=plt.figure(); ax=f.add_subplot(111)
    ax.plot(time, taunew, 'k-', label='new taup')
    ax.plot([min(time), max(time)], [taup*1e-3, taup*1e-3], 'r-', label='TAUPH set')
    ax.legend(loc='best')
   

    uf_d = {'pre': 'OMF',
            'ext': 'TPI',
            'shot': 58823,
            'grid': {'X': {'lbl': 'Time'+' '*19+'SEC'+' '*4,
                           'arr': time}},
                           'data': {'lbl': 'Confinement_time       SEC',
                                    'arr': taunew}}
                
    uf.WU(uf_d, udir='/home/vallar')
