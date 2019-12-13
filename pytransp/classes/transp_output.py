"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
from __future__ import print_function

import netCDF4 as nc
import numpy as np
import pytransp.trutils.transp_utils as tu
import pytransp.plot.plot_input as plot_input
import pytransp.plot.plot_eq as plot_eq

class transp_output:
    """Class for 1d output
       
    Parameters:
        None
    Attributes:
        None
    Note:

    """
    def __init__(self, fname):
        """Initialization from netCDF file
                
        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self.fname = fname
        self.file = nc.Dataset(self.fname)
        self.runid = self.fname[-12:-4]
        self.run = self.runid[0:5]
        self._set_vars()
        self._calculate_all()
        self.flag_beam=0
        
    def _set_vars(self):
        """Initialization
        
        Hidden function to initialize variables

        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self._global_vars()
        self._imp_vars()
        self._kinetic_vars()

    def _calculate_all(self):
        """first computations
        
        Hidden function to initialize and calculate variables
        i.e. vol-averaged temperatures and densities, ecc.

        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self._average_kin()
        
    def _global_vars(self):
        """Global values

        Gathers global values, i.e. time and nt, rho grid, shot, run, tokamak
        usually the 2D arrays are stored as (time, rho)
        
        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self.t  = self.file['TIME'][:]
        self.nt = len(self.t)
        self.darea = np.array(self.file.variables['DAREA'][:])*1e-4 # in m2
        self.dvol  = np.array(self.file.variables['DVOL'][:])*1e-6 #in m3
        self.rho = np.array(self.file.variables['X'][:])
        
    def _imp_vars(self):
        """Impurities

        Gathers info on impurities

        Parameters:
            None
        Attributes:
            None
        Note:

        """
        keys     = ['Zeff', 'A', 'Z', 'nimp']
        varnames = ['ZEFFP', 'AIMP', 'ZIMPS_TOK', 'NIMP']
        self.imp_names, self.imp_vars = tu._fill_dict(self.file, keys, varnames)        
        
    def _kinetic_vars(self):
        """Kinetic variables

        Gathers kinetic variables (ne,ni,nimp,te,ti) from netcdf file

        Parameters:
            None
        Attributes:
            None
        Note:

        """
        keys     = ['ne', 'ni', 'te', 'ti']
        varnames = ['NE', 'NI', 'TE', 'TI']
        self.kin_names = dict.fromkeys(keys)
        self.kin_vars  = dict.fromkeys(keys)
        self.kin_names, self.kin_vars = \
            tu._fill_dict(self.file, keys, varnames)
        self.kin_vars['ne'] *= 1e6
        self.kin_vars['ni'] *= 1e6


    def _average_kin(self):
        """ Average of kin. vars

        Calculates the vol-average of densities and temperatures

        Parameters:
            None
        Attributes:
            None    
        Note:

        """
        self.ne_vavg = np.zeros(self.nt, dtype=float)
        self.ni_vavg = np.zeros(self.nt, dtype=float)
        self.Te_vavg = np.zeros(self.nt, dtype=float)
        self.Ti_vavg = np.zeros(self.nt, dtype=float)
        self.ip      = np.zeros(self.nt, dtype=float)
        for i in range(self.nt):
            dvol = self.dvol[i,:]; darea=self.darea[i,:]
            self.ne_vavg[i] = np.dot(self.kin_vars['ne'][i,:],dvol)/np.sum(dvol)
            self.ni_vavg[i] = np.dot(self.kin_vars['ni'][i,:],dvol)/np.sum(dvol)
            self.Te_vavg[i] = np.dot(self.kin_vars['te'][i,:],dvol)/np.sum(dvol)
            self.Ti_vavg[i] = np.dot(self.kin_vars['ti'][i,:],dvol)/np.sum(dvol)
            self.ip[i]      = np.dot(self.file.variables['CUR'][i,:]*1e4, darea)
        self.ne_mean = np.mean(self.ne_vavg)
        self.ni_mean = np.mean(self.ni_vavg)
        self.Te_mean = np.mean(self.Te_vavg)
        self.Ti_mean = np.mean(self.Ti_vavg)


    def plot_input(self, time=[0]):
        """
        """
        plot_input.plot_input(self, time)
        
    def plot_input_1d(self, time=[0], axne=0, axTe=0, ls='-'):
        """
        """
        plot_input.plot_input_1d(self, time, axne, axTe, ls)
        
    def plot_input_prof(self, time=[0]):
        """
        """
        plot_input.plot_input_prof(self, time)


    def _gather_equilibrium(self):
        """get data for eq
        
        Stores data coming from equilibrium, i.e.
        j, F, p, q, psi, phi
        
        """
        keys     = ['j', 'joh', 'jbs', 'p', 'pth', 'pnth', 'pol_flux', 'tor_flux', 'q', 'f']
        varnames = ['CUR', 'CUROH', 'CURBS', 'PMHD_IN', 'PMHDT_IN', 'PMHDF_IN', 'PLFLX', 'TRFLX','Q', 'BZXR']
        self.eq_names = dict.fromkeys(keys)
        self.eq_vars  = dict.fromkeys(keys)
        self.eq_names, self.eq_vars = \
            tu._fill_dict(self.file, keys, varnames)
        self.eq_vars['j']*=1e4
        self.eq_vars['joh']*=1e4
        self.eq_vars['jbs']*=1e4
        self.eq_vars['f']*=1e-2
        
    def plot_equilibrium(self, time=[0],f=0):
        """
        """
        try:
            self.eq_vars['j'].mean()
        except:
            self._gather_equilibrium()
        plot_eq.plot_eq(self, time,f)
