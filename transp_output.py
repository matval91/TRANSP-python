#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

col = ['k','r','b','m','g','c']

class transp_output:
    """
    """
    def __init__(self, fname):
        """
        """
        
        self.fname = fname
        self.file = nc.Dataset(self.fname)
        self._global_vars()
        self._kinetic_vars()
        self._nb_vars()
        
    def _global_vars(self):
        """
        Gathers global values, i.e. time and nt, rho grid, shot, run, tokamak
        usually the 2D arrays are stored as (time, rho)
        """
        self.t  = self.file['TIME'][:]
        self.nt = len(self.t)
        self.darea = self.file.variables['DAREA'][:]
        self.dvol  = self.file.variables['DVOL'][:]
        self.rho = np.linspace(0, 1, num=len(self.darea[0,:]), dtype=float)

        self._imp_vars()
        
    def _imp_vars(self):
        """
        Gathers info on impurities
        """
        keys     = ['Zeff', 'A', 'Z', 'nimp']
        varnames = ['ZEFFP', 'AIMP', 'ZIMPS_TOK', 'NIMP']
        self.imp_names = dict.fromkeys(keys)
        self.imp_vars  = dict.fromkeys(keys)
        self.imp_names, self.imp_vars = \
            self._fill_dict(keys, varnames, self.imp_names, self.imp_vars)        
        
    def _kinetic_vars(self):
        """
        Gathers kinetic variables (ne,ni,nimp,te,ti) from netcdf file
        """
        keys     = ['ne', 'ni', 'te', 'ti']
        varnames = ['NE', 'NI', 'TE', 'TI']
        self.kin_names = dict.fromkeys(keys)
        self.kin_vars  = dict.fromkeys(keys)
        self.kin_names, self.kin_vars = \
            self._fill_dict(keys, varnames, self.kin_names, self.kin_vars)

    def _nb_vars(self):
        self._nb_in_vars()
        self._nb_out_vars()

    def _nb_in_vars(self):
        """
        Gathers data in input to TRANSP: 
        beam energy, beam power, beam species, beam geometry
        """
        keys     = ['P', 'E']
        varnames = ['PBINJ_D', 'EINJAV_D' ]
        self.nb_in_names = dict.fromkeys(keys)
        self.nb_in_vars  = dict.fromkeys(keys)
        self.nb_in_names, self.nb_in_vars = \
            self._fill_dict(keys, varnames, self.nb_in_names, self.nb_in_vars)        

    def _nb_out_vars(self):
        """
        Gathers NB variables (nfast, pfast, pe, pi, nbcd) from netcdf file
        """
        keys     = ['st']
        varnames = ['PBSHINE_D']
        self.nb_ioniz_names = dict.fromkeys(keys)
        self.nb_ioniz_vars  = dict.fromkeys(keys)
        self.nb_ioniz_names, self.nb_ioniz_vars = \
            self._fill_dict(keys, varnames, self.nb_ioniz_names, self.nb_ioniz_vars)         
        
        keys     = ['orbloss', 'cx1', 'cx2', 'pi', 'pe',\
                    'th', 'nbcd']
        varnames = ['BPLIM_D', 'BPCXX_D', 'BPCXI_D', 'PBI_D', 'PBE_D', \
                    'PBTH_D', 'CURB']
        self.nb_FKP_names = dict.fromkeys(keys)
        self.nb_FKP_vars  = dict.fromkeys(keys)
        self.nb_FKP_names, self.nb_FKP_vars = \
            self._fill_dict(keys, varnames, self.nb_FKP_names, self.nb_FKP_vars) 
            
        
    def _fill_dict(self, keys, varnames, name_dict, var_dict):
        """
        Fills dictionaries
        """
        name_dict = dict.fromkeys(keys)
        var_dict  = dict.fromkeys(keys)
        
        for ii,kk in enumerate(keys):
            name_dict[kk] = varnames[ii]
            tmpdata = self.file.variables[name_dict[kk]][:]
            var_dict[kk] = tmpdata
            
        return name_dict, var_dict
    
    
    def compute_Ec(self):
        """
        Ec = 14.8*Te*(A^{3/2}/ne*sum(nj*Zj/Aj))^2/3
        """
        try:
            self.imp_vars.keys()
        except:
            self._imp_vars()
        print "BEAM SPECIES SET TO D"
        A=2
        Aj = [2, self.imp_vars['A']]
        Zj = [1, self.imp_vars['Z']]
        nj = [self.kin_vars['ni'], self.imp_vars['nimp']]
        ne = self.kin_vars['ne']; Te = self.kin_vars['te']
        
        sum_facts = [nj[i]*Zj[i]**2/Aj[i] for i in range(len(Aj))]
        summ = np.sum(np.asarray(sum_facts), axis=0)
        second_factor = (A**1.5/ne*summ)**0.66
        Ec_prof = 14.8*Te*second_factor
        self.Ec = Ec_prof
        
    def compute_tau(self):
        """
        # A plasma
        tau_s=ts/3*ln(1+(E/EC)**1.5)
        ts=6.27e8*A Te[eV]**1.5/(Z**2 ne(cm^-3) lnLambda)
        valid for Z=1
        """
        try:
            self.Ec.mean()
        except:
            self.compute_Ec()

        # Checks density normalisation (it needs m^-3, not cm^-3)

        A = [2, self.imp_vars['A']]
        ne = self.kin_vars['ne']; Te = self.kin_vars['te']
        Ec = self.Ec; E0 = self.nb_in_vars['E']
        lnlambda=17.
        if np.mean(ne)/1e18 >1:
            ne = ne*1e-6
        else:
            ne=ne
        
        ts=6.27e8*A[0]*np.divide(np.power(Te,1.5),(1.*lnlambda*ne))
        E0=np.transpose(np.tile(E0, (np.shape(Ec)[1], 1)))
        taus=ts/3.*(1+(E0/Ec)**1.5)
        self.taus = taus

        
    def compute_Gi(self):
        """
        Gi = Ec/E0 int[0, E0/Ec]dy/(1+y**1.5)
        """
        try:
            self.Ec.mean()
        except:
            self.compute_Ec()
        Ec = self.Ec; E0 = self.nb_in_vars['E']
        gi = np.zeros(np.shape(Ec), dtype=float)
        for time_ind in range(np.shape(Ec)[0]):
            timeslice = Ec[time_ind,:]
            if E0[time_ind]==0:
                gi[time_ind,:] = np.zeros(np.shape(Ec)[1])
            else:
                for rho_ind in range(np.shape(Ec)[1]):
                    xarray = np.linspace(0,E0[time_ind]/timeslice[rho_ind], num=100, dtype=float)
                    gi[time_ind,rho_ind] = timeslice[rho_ind]/E0[time_ind]*np.trapz(1/(1+xarray**1.5), xarray)
        self.gi = gi
        
    def plot_deposition(self):
        """
        Plots deposition to ions, electrons, current drive
        """
        try:
            self.nb_FKP_vars['pi'].mean()
        except:
            print "No Pi data"
            return
        if self.nt > 5:
            ind = np.trunc(np.linspace(0, self.nt-1, 5))
        ind = ind.astype(int)
        print ind
        f = plt.figure()
        axe = f.add_subplot(121)
        axi = f.add_subplot(122)
        for i in ind:
            axi.plot(self.rho, self.nb_FKP_vars['pi'][i,:], col[i], label=r't = '+str(self.t[i]))
            axe.plot(self.rho, self.nb_FKP_vars['pe'][i,:], col[i], label=r't = '+str(self.t[i]))
        axe.set_xlabel(r'$\rho$'); axe.set_ylabel(r'Power deposition [$W/m^3$]')
        axi.set_xlabel(r'$\rho$'); axi.set_ylabel(r'Power deposition [$W/m^3$]')

        axe.legend(loc='best')
        f.tight_layout()
        plt.show()