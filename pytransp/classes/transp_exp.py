#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

from pytransp.classes.transp_heating import transp_heating
import pytransp.trutils.transp_utils as tu
import utils.plot_utils as au

class transp_exp(transp_heating):
    """Class for 1d output
       
    Parameters:
        None
    Attributes:
        None
    Note:

    """
    def __init__(self, fname):
        """Initialization from transp_output class
                
        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self.fname = fname
        transp_heating.__init__(self,fname)
        self._calculate_n0()
        self._calculate_wmhd()
        self._performance_vars()

    def _performance_vars(self):
        """
        Gathers performance variables (jbootstrap, beta, li) from netcdf file
        ), (u'CURBS', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBS(TIME3, X)
            units: AMPS/CM2
            long_name: BOOTSTRAP CURRENT
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSNEO', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSNEO(TIME3, X)
            units: AMPS/CM2
            long_name: NEO-gk Bootstrap Current
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSWNC', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSWNC(TIME3, X)
            units: AMPS/CM2
            long_name: NCLASS Bootstrap Current
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSEPS', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSEPS(TIME3, X)
            units: AMPS/CM2
            long_name: Aspect Ratio Bootstrap Current
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSSAU', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSSAU(TIME3, X)
            units: AMPS/CM2
            long_name: Sauter Bootstrap Current as Used
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSSAU0', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSSAU0(TIME3, X)
            units: AMPS/CM2
            long_name: Sauter Bootstrap Current Original Form
        unlimited dimensions:
        current shape = (101, 40)
        filling off
        ), (u'CURBSSAU1', <type 'netCDF4._netCDF4.Variable'>
        float32 CURBSSAU1(TIME3, X)
            units: AMPS/CM2
            long_name: Sauter Bootstrap Current CS Chang Form
        unlimited dimensions:
        current shape = (101, 40)
        """        
        keys     = ['jboot_no', 'jboot', 'jbootneogk', 'jbootneo', 'jbootsauor']
        varnames = ['CURBS', 'CURBSSAU', 'CURBSNEO', 'CURBSWNC', 'CURBSSAU0']
        self.perf_names = dict.fromkeys(keys)
        self.perf_vars  = dict.fromkeys(keys)
        self.perf_names, self.perf_vars = tu._fill_dict(self.file,keys, varnames)
        for el in self.perf_vars:
            self.perf_vars[el] *= 1e4
        self.jboot=np.copy(self.t)
        
        for i,t in enumerate(self.t):
            self.jboot[i]=np.dot(self.darea[i,:], self.perf_vars['jboot'][i,:])

    def n0_core(self):
        """ plot n0 at core

        plot n0 at core
       
        Parameters:
            None
        Attributes:
            None    
        Note:
        """   
        v_source = self.file.variables['DN0VD'][:]
        w_source = self.file.variables['DN0WD'][:]

        CXfastn = self.file.variables['N0BCXD0'][:]
        first_fastn = self.file.variables['N0BD0'][:]
        halob = self.file.variables['N0BH_D'][:]
        Drecy = self.file.variables['N0RC_D_D'][:]
        Dflow= self.file.variables['N0GF_D_D'][:]
        Dsflow= self.file.variables['N0SGF_D'][:]
        Dnrecy= self.file.variables['N0SRC_D'][:]
        Drec= self.file.variables['N0V0_D'][:]
         
        tot_source = v_source+w_source+first_fastn+CXfastn+Drecy+Dflow+halob+Dsflow+Dnrecy+Drec

        f=plt.figure(); ax=f.add_subplot(111)
        ax.plot(self.t, tot_source[:, 0]*1e6, 'k', lw=2.3, label='tot')
        ax.plot(self.t, v_source[:,0]*1e6, 'b', lw=2.3, label='Volume')
        ax.plot(self.t, halob[:,0]*1e6, 'm', lw=2.3, label='halo')
        ax.plot(self.t, first_fastn[:,0]*1e6, 'c', lw=2.3, label='fastn')
        ax.plot(self.t, w_source[:,0]*1e6, 'r', lw=2.3, label='wall')
        ax.plot(self.t, Drecy[:,0]*1e6+Dnrecy[:,0]*1e6+Drec[:,0]*1e6, 'g', lw=2.3, label='Recy')
        ax.plot(self.t, CXfastn[:,0]*1e6,'y', lw=2.3, label='CX')
        ax.set_xlabel(r't [s]'); ax.set_ylabel(r'n0 [1/m3]')
        ax.legend(loc='best'); ax.grid('on')
        plt.show()

    def plot_n0_edge(self):
        """ plot n0 at lcfs

        plot n0 at lcfs
       
        Parameters:
            None
        Attributes:
            None    
        Note:
        """  
        self._calculate_n0(plot=1)


    def _calculate_n0(self, plot=0):
        """ plot n0 at lcfs
        
        plot n0 at lcfs
        
        Parameters:
            None
        Attributes:
            None    
        Note:
        """   
        v_source = self.file.variables['DN0VD'][:]
        w_source = self.file.variables['DN0WD'][:]
#        Drecy = self.file.variables['N0RC_D_D'][:]
#        Dflow= self.file.variables['N0GF_D_D'][:]
#        Dsflow= self.file.variables['N0SGF_D'][:]
#        Dnrecy= self.file.variables['N0SRC_D'][:]
#        Drec= self.file.variables['N0V0_D'][:]
#        
        CXfastn = self.file.variables['N0BCXD0'][:]
        first_fastn = self.file.variables['N0BD0'][:]
#        halob = self.file.variables['N0BH_D'][:]
        n0fast = first_fastn+CXfastn
        tot_source = v_source+w_source#+Drecy+Dflow+Dsflow+Dnrecy+Drec
        tot_source += n0fast
        
        self.n0_tot = tot_source*1e6 #in m^-3
        if plot==1:
            au.common_style()
            f=plt.figure(); ax=f.add_subplot(111)
            ax.plot(self.t, tot_source[:, -1]*1e6, 'k', lw=2.3, label='tot')
            ax.plot(self.t, v_source[:,-1]*1e6, 'b', lw=2.3, label='Volume')
            ax.plot(self.t, Drecy[:,-1]*1e6+Dnrecy[:,-1]*1e6+Drec[:,-1]*1e6, 'g', lw=2.3, label='Recy')
            ax.plot(self.t, w_source[:,-1]*1e6, 'r', lw=2.3, label='wall')
            ax.plot(self.t, n0fast[:,-1]*1e6, 'm', lw=2.3, label='fast')
            
            
#            ax.plot(self.t, halob[:,-1]*1e6, 'm', lw=2.3, label='halo')
#            ax.plot(self.t, first_fastn[:,-1]*1e6, 'c', lw=2.3, label='fastn')
#            ax.plot(self.t, CXfastn[:,-1]*1e6,'y', lw=2.3, label='CX')
            ax.set_xlabel(r't [s]'); ax.set_ylabel(r'n0 [1/m3]')
            ax.legend(loc='best'); ax.grid('on')
            plt.show()

    def _calculate_wmhd(self):
        """
        Calculate the perpendicular stored energy
        wi, we, wimp = perpendicular energies of bulk ions, electrons and impurities
        wth= only thermal component of the energy
        wprp = perpendicular component of fast ion energy
        """
        if self.flag_beam==1:
            wi  = np.multiply(self.kin_vars['ni']-self.nb_FKP_vars['n'],self.kin_vars['ti'])*1.602e-19*1.5
            #wi  = np.multiply(self.file.variables['ND'][:,:]*1e6-self.nb_FKP_vars['n'],self.kin_vars['ti'])*1.602e-19*1.5
        else: 
            wi  = np.multiply(self.kin_vars['ni'],self.kin_vars['ti'])*1.602e-19*1.5
            #wi  = np.multiply(self.file.variables['ND'][:,:],self.kin_vars['ti'])*1.602e-19*1.5

        we  = np.multiply(self.kin_vars['ne'],self.kin_vars['te'])*1.602e-19*1.5
        wim = np.multiply(self.imp_vars['nimp'],self.kin_vars['ti'])*1.602e-19*1.5
        wth = we+wi+wim
        wtot_dens = self.file.variables['UTOTL'][:]
        
        self.wth=np.zeros(len(self.t)); self.wtot=np.copy(self.wth); self.wprp=np.copy(self.wth); self.we=np.copy(self.wth); self.wi=np.copy(self.wth)
        for i in range(len(self.t)):
            self.wth[i]  = np.dot(wth[i], self.dvol[i,:])
            self.we[i] = np.dot(we[i], self.dvol[i,:])
            self.wi[i] = np.dot(wi[i], self.dvol[i,:])
            if self.flag_beam==1:
                self.wprp[i] = np.dot(self.file.variables['UBPRP'][i]*1e6, self.dvol[i,:])
            self.wtot[i] = np.dot(wtot_dens[i,:]*1e6, self.dvol[i,:])

        #self.wtot=self.wth+1.5*self.wprp
