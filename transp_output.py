#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, interactive
import math, collections
import scipy.interpolate as interp
import glob, os, shutil

col = ['k','r','b','m','g','c']

cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)

class transp_output_1d:
    """
    """
    def __init__(self, fname):
        """
        """
        
        self.fname = fname
        self.file = nc.Dataset(self.fname)
        self._set_vars()
        
    def _set_vars(self):
        """
        Gathers all the values needed
        """
        self._global_vars()
        self._kinetic_vars()
        self._nb_vars() 
        self._performance_vars()
        
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
        self.kin_vars['ne'] *= 1e6
        self.kin_vars['ni'] *= 1e6

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
        
        keys     = ['n','orbloss', 'cx1', 'cx2', 'pi', 'pe',\
                    'th', 'nbcd', 'cxprof', 'pMC']
        varnames = ['BDENS_D','BPLIM_D', 'BPCXX_D', 'BPCXI_D', 'PBI_D', 'PBE_D', \
                    'PBTH_D', 'CURB', 'PBCX_D', 'PBDEPMC_D']
        self.nb_FKP_names = dict.fromkeys(keys)
        self.nb_FKP_vars  = dict.fromkeys(keys)
        self.nb_FKP_names, self.nb_FKP_vars = \
            self._fill_dict(keys, varnames, self.nb_FKP_names, self.nb_FKP_vars) 
        self.nb_FKP_vars['n'] *= 1e6
     
    def _performance_vars(self):
        """
        Gathers performance variables (jbootstrap, beta, li) from netcdf file
        """        
        keys     = ['jboot']
        varnames = ['CURBS']
        self.perf_names = dict.fromkeys(keys)
        self.perf_vars  = dict.fromkeys(keys)
        self.perf_names, self.perf_vars = \
            self._fill_dict(keys, varnames, self.perf_names, self.perf_vars)
            
            
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
    
    
    def calculate_Ec(self):
        """
        Ec = 14.8*Te*(A^{3/2}/ne*sum(nj*Zj/Aj))^2/3
        """
        try:
            self.imp_vars.keys()
        except:
            self._imp_vars()
        print "BEAM SPECIES SET TO D"
        A=2
        Aj = [2, self.imp_vars['A'].mean()]
        Zj = [1, self.imp_vars['Z'].mean()]
        nj = [self.kin_vars['ni'], self.imp_vars['nimp']]
        ne = self.kin_vars['ne']; Te = self.kin_vars['te']
        self.Ec = np.zeros((self.nt, len(self.rho)) , dtype=float)
        self.ec_mean = np.zeros(self.nt, dtype=float)
        
        for it in range(self.nt):
            sum_facts = [nj[i][it,:]*(Zj[i]**2.)/Aj[i] for i in range(len(Aj))]
            summ = np.sum(np.asarray(sum_facts), axis=0)
            second_factor = (A**1.5/ne[it,:]*summ)**0.66
            Ec_prof = 14.8*Te[it,:]*second_factor
            self.Ec[it,:] = Ec_prof
            self.ec_mean[it] = np.trapz(Ec_prof, self.rho)
            
    def calculate_tau(self):
        """
        # A plasma
        tau_s=ts/3*ln(1+(E/EC)**1.5)
        ts=6.27e8*A Te[eV]**1.5/(Z**2 ne(cm^-3) lnLambda)
        valid for Z=1
        """
        try:
            self.Ec.mean()
        except:
            self.calculate_Ec()

        # Checks density normalisation (it needs m^-3, not cm^-3)

        A = [2, self.imp_vars['A'].mean()]
        ne = self.kin_vars['ne']; Te = self.kin_vars['te']
        Ec = self.Ec; E0_arr = self.nb_in_vars['E']
        lnlambda=17.
        self.taus = np.zeros((self.nt, len(self.rho)) , dtype=float)
        self.taus_mean = np.zeros(self.nt, dtype=float)
        
        if np.mean(ne)/1e18 >1:
            ne = ne*1e-6
        else:
            ne=ne
            
        for it in range(self.nt):
            ts=6.27e8*A[0]*np.divide(np.power(Te[it,:],1.5),(1.*lnlambda*ne[it,:]))
            E0=np.transpose(np.tile(E0_arr[it], (len(self.rho), 1)))
            taus=ts/3.*(1+(E0/Ec[it,:])**1.5)
            self.taus[it,:] = taus
            self.taus_mean[it] = np.trapz(taus, self.rho)

        
    def calculate_gi(self):
        """
        Gi = Ec/E0 int[0, E0/Ec]dy/(1+y**1.5)
        """
        try:
            self.Ec.mean()
        except:
            self.calculate_Ec()
        Ec = self.Ec; E0 = self.nb_in_vars['E']
        self.gi = np.zeros(np.shape(Ec), dtype=float)
        self.gi_mean = np.zeros(self.nt, dtype=float)
        for time_ind in range(np.shape(Ec)[0]):
            if E0[time_ind] == 0.:
                continue
            timeslice = Ec[time_ind,:]
            for rho_ind in range(np.shape(Ec)[1]):
                if timeslice[rho_ind] == 0:
                    self.gi[time_ind, rho_ind] = 0.
                else:
                    x_array = np.linspace(0,E0[time_ind]/timeslice[rho_ind], num=100, dtype=float)
                    self.gi[time_ind,rho_ind] = timeslice[rho_ind]/E0[time_ind]*np.trapz(1/(1+x_array**1.5), x_array)
            self.gi_mean[time_ind] = np.trapz(self.gi[time_ind,:], self.rho)   
        
        
    def plot_input(self):
        """
        Plots input quantities, such as ne, ni, Te, Ti
        """
        try:
            self.kin_vars['ne'].mean()
        except:
            print "No ne data"
            return
        if self.nt > 5:
            ind = np.trunc(np.linspace(0, self.nt-1, 5))
        ind = ind.astype(int)
        
        f = plt.figure()
        axne = f.add_subplot(221)
        axni = f.add_subplot(222, sharey=axne)
        axTe = f.add_subplot(223)
        axTi = f.add_subplot(224, sharey=axTe)
        
        for i, el in enumerate(ind):
            axne.plot(self.rho, self.kin_vars['ne'][el,:], col[i], label=r't = '+str(self.t[el]))
            axni.plot(self.rho, self.kin_vars['ni'][el,:], col[i], label=r't = '+str(self.t[el]))
            axTe.plot(self.rho, self.kin_vars['te'][el,:]*1e-3, col[i], label=r't = '+str(self.t[el]))
            axTi.plot(self.rho, self.kin_vars['ti'][el,:]*1e-3, col[i], label=r't = '+str(self.t[el]))

        axne.set_xlabel(r'$\rho$'); axne.set_ylabel(r'$n_e$ [$1/m^3$]')
        axTe.set_xlabel(r'$\rho$'); axTe.set_ylabel(r'$T_e$ [$keV$]')

        axni.set_xlabel(r'$\rho$'); axni.set_ylabel(r'$n_i$ [$1/m^3$]')
        axTi.set_xlabel(r'$\rho$'); axTi.set_ylabel(r'$T_i$ [$keV$]')

        axne.legend(loc='best'), axTe.legend(loc='best')
        f.suptitle(self.fname)
        f.tight_layout()
        plt.show()    
    
    def plot_deposition(self):
        """
        Plots deposition to ions, electrons
        """
        try:
            self.nb_FKP_vars['pi'].mean()
        except:
            print "No Pi data"
            return
        if self.nt > 5:
            ind = np.trunc(np.linspace(0, self.nt-1, 5))
        ind = ind.astype(int)
        f = plt.figure()
        axe = f.add_subplot(121)
        axi = f.add_subplot(122)
        for i, el in enumerate(ind):
            axi.plot(self.rho, self.nb_FKP_vars['pi'][el,:], col[i], label=r't = '+str(self.t[el]))
            axe.plot(self.rho, self.nb_FKP_vars['pe'][el,:], col[i], label=r't = '+str(self.t[el]))
        axe.set_xlabel(r'$\rho$'); axe.set_ylabel(r'$P_e$ [$W/m^3$]')
        axi.set_xlabel(r'$\rho$'); axi.set_ylabel(r'$P_i$ [$W/m^3$]')

        axe.legend(loc='best')
        f.suptitle(self.fname)
        f.tight_layout()
        plt.show()
        
        
    def plot_NBCD(self):
        """
        Plot current drive
        """
        try:
            self.nb_FKP_vars['nbcd'].mean()
        except:
            print "No nbcd data"
            return
        if self.nt > 5:
            ind = np.trunc(np.linspace(0, self.nt-1, num=5))
        ind = ind.astype(int)
        self.calculate_scalars()

        f = plt.figure()
        axj  = f.add_subplot(121)
        axcd = f.add_subplot(122)
        axcd.plot(self.t, self.nbcd*1e-3, 'k', lw=3.)
        for i, el in enumerate(ind):
            axj.plot(self.rho, self.nb_FKP_vars['nbcd'][el,:], col[i], label=r't = '+str(self.t[el]))
        
        axj.set_xlabel(r'$\rho$'); axj.set_ylabel(r'j shielded [$A/m^2$]')
        axcd.set_xlabel(r't [s]'); axcd.set_ylabel(r'Driven current [kA]')
        axj.legend(loc='best')
        f.tight_layout(); f.suptitle(self.fname)

        plt.show()
     
     
    def average_kin(self):
        """
        Calculates the average of densities and temperatures
        """
        self.ne_mean = np.zeros(self.nt, dtype=float)
        self.ni_mean = np.zeros(self.nt, dtype=float)
        self.Te_mean = np.zeros(self.nt, dtype=float)
        self.Ti_mean = np.zeros(self.nt, dtype=float)
        self.nf_mean = np.zeros(self.nt, dtype=float)
       
        for i in range(self.nt):
            self.ne_mean[i] = np.trapz(self.kin_vars['ne'][i,:], self.rho)
            self.ni_mean[i] = np.trapz(self.kin_vars['ne'][i,:], self.rho)
            self.Te_mean[i] = np.trapz(self.kin_vars['ne'][i,:], self.rho)
            self.Ti_mean[i] = np.trapz(self.kin_vars['ne'][i,:], self.rho)
            self.nf_mean[i] = np.trapz(self.nb_FKP_vars['n'][i,:], self.rho)
        
    def _calculate_scalarpower(self):
        """
        Calculates scalar variables
        """
        #CD
        self.pe = np.zeros(len(self.t), dtype=float)
        self.pi = np.zeros(len(self.t), dtype=float)
        self.pth = np.zeros(len(self.t), dtype=float)
        self.pcxprof = np.zeros(len(self.t), dtype=float)
        for i in range(len(self.t)):
            self.pe[i] = np.dot(self.nb_FKP_vars['pe'][i,:], self.dvol[i,:])
            self.pi[i] = np.dot(self.nb_FKP_vars['pi'][i,:], self.dvol[i,:])
            self.pth[i] = np.dot(self.nb_FKP_vars['th'][i,:], self.dvol[i,:])
            self.pcxprof[i] = np.dot(self.nb_FKP_vars['cxprof'][i,:], self.dvol[i,:])
        self.pcx = self.nb_FKP_vars['cx1']+self.nb_FKP_vars['cx2']
        self.pol = self.nb_FKP_vars['orbloss']
        
        self.psh = self.nb_ioniz_vars['st']
            
    def power_balance(self):
        """
        """            
        self._calculate_scalarpower()
        pbal = np.mean(self.nb_in_vars['P']-self.psh - \
            self.pe-self.pi-self.pth-self.pol-self.pcx)*1e-6
        print "POWER BALANCE:", pbal
        print "Pbalance/Pin ", pbal/np.mean(self.nb_in_vars['P'])*1e6
        print "psh+pol ", np.mean(self.pol+self.psh)*1e-6
        print "pe+pi+pth ", np.mean(self.pe+self.pi+self.pth)*1e-6
        print "pcx " , np.mean(self.pcxprof)*1e-6
        
    def calculate_scalars(self):
        """
        """
        self.nbcd = np.zeros(len(self.t), dtype=float)
        self.jboot = np.zeros(len(self.t), dtype=float)
        for i in range(len(self.t)):
            self.nbcd[i] = np.dot(self.nb_FKP_vars['nbcd'][i,:], self.darea[i,:])
            self.jboot[i] = np.dot(self.perf_vars['jboot'][i,:], self.darea[i,:])
        self._calculate_scalarpower()
        self.average_kin()
        
    def plot_performance(self):
        """
        """
        try:
            self.perf_vars['jboot'].mean()
        except:
            print "No bootstrap data"
            return
        self.calculate_scalars()
        f = plt.figure()
        axcd = f.add_subplot(121)
        axcd.plot(self.t, self.jboot*1e-3, 'k', lw=3.)
        axcd.set_xlabel(r't [s]'); axcd.set_ylabel(r'Bootstrap current [kA]')
        f.tight_layout(), f.suptitle(self.fname)
        plt.show()


    def calc_eff(self):                
        """
        Method to calculate the current-drive efficiency from the beam:
        eta = R0*n_e*I_CD/P
        """                
        try:
            self.nbcd.mean()
        except:
            self.calculate_scalars()
        #Computing the power coupled to plasma: Pini-Pend+Pres
        P = self.nb_in_vars['P']
        #ind = P>0.
        R0 = self.file.variables['RAXIS'][:]*0.01
        # average density in 10^20
        ne_avg = self.ne_mean[:]*1e-20
        #############################
        self.eff = R0*ne_avg*np.abs(self.nbcd)/P
        self.eff[np.isinf(self.eff)] = 0.
        self.eff[np.isnan(self.eff)] = 0.
#        print "R0 ", R0, 'm'
#        print "ne avg", ne_avg, '10e20 m^-3'
#        print "Ip ", self.nbcd, 'A'
#        print "P  ", P, 'W'
        print "CD efficiency: ", self.eff, " [10^20 A / (W m^2)]"


class transp_deposition:
    """
    Class to analyse the denposition profiles from NUBEAM
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_r_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_r_D_MCBEAM: R at deposition, cm
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_z_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_z_D_MCBEAM: Z at deposition, cm
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_rgc_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_rgc_D_MCBEAM: R at GC, cm
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_zgc_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_zgc_D_MCBEAM: Z at GC, cm
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_xksid_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_xksid_D_MCBEAM: pitch angle, V_II/V
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_einj_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_einj_D_MCBEAM: Energy, eV
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_wght_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_wght_D_MCBEAM: weight deposited ptcl
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_zeta_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_zeta_D_MCBEAM: toroidal angle at deposition, deg
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    float64 bs_time_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_time_D_MCBEAM: TIME at deposition, sec
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    
    <type 'netCDF4._netCDF4.Variable'>
    int32 bs_ib_D_MCBEAM(nbsdep_D_MCBEAM)
        bs_ib_D_MCBEAM: BEAM ID
    unlimited dimensions: 
    current shape = (50478,)
    filling off
    """
    def __init__(self, inf_name):
        """
        Gathers variables
        """
        self.infile_n = inf_name
        self.file = nc.Dataset(inf_name)
        
        keys = ['R','z', 'Rgc','zgc', 'pitch','E','weight','phi','time', 'beamid']
        varnames = ['bs_r_D_MCBEAM', 'bs_z_D_MCBEAM', 'bs_rgc_D_MCBEAM',\
            'bs_zgc_D_MCBEAM', 'bs_xksid_D_MCBEAM', 'bs_einj_D_MCBEAM', \
            'bs_wght_D_MCBEAM', 'bs_zeta_D_MCBEAM', 'bs_time_D_MCBEAM', \
            'bs_ib_D_MCBEAM']
        self.data_i_names = dict.fromkeys(keys)
        self.data_i  = dict.fromkeys(keys)

        self.data_i_names, self.data_i = self._fill_dict(keys, varnames, \
                self.data_i_names, self.data_i) 
        self.npart = len(self.data_i['R'])
        
        #Convert phi from deg to rad
        self.data_i['phi']=self.data_i['phi']*math.pi/180.

        #convert the distances from cm to m
        for i in ['R','z', 'Rgc', 'zgc']:
            self.data_i[i] *= 0.01
            


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



    def plot_RZpart(self):
        """
        Method to plot R vs z of the ionised particles, useful mostly with bbnbi
        """
        x=self.data_i['R']
        y=self.data_i['z']

        f=plt.figure()
        ax = f.add_subplot(111)
        #hb = ax.hexbin(x, y, gridsize=100, cmap=my_cmap)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap)
        cb = f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('R')
        ax.set_ylabel('z')

        if len(self.R_w)==0:
            self._readwall()
        ax.plot(self.R_w , self.z_w, 'k', linewidth=3)

        ax.axis('equal')
        ax.axis([min(self.R_w), max(self.R_w), min(self.z_w), max(self.z_w)])
        #ax.set_xrange([min(self.R_w), max(self.R_w)])
        #ax.set_yrange([min(self.z_w), max(self.z_w)])                   
        plt.show()
        
    def plot_XYpart(self):
        """
        Method to plot XY of ionisation, without difference between the beams
        """

        x=np.zeros(self.npart)
        y=np.zeros(self.npart)
        
        for i, el in enumerate(self.data_i['R']):
            x[i]=el*math.cos(self.data_i['phi'][i])
            y[i]=el*math.sin(self.data_i['phi'][i])
            
        f=plt.figure()
        ax=f.add_subplot(111)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap, normed=True)
        cb = f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        theta=np.arange(0,6.3,0.02*6.28)
        #if len(self.R_w)==0:
        #    ax.plot(np.min(self.R_w)*np.cos(theta) , np.min(self.R_w)*np.sin(theta), 'k', linewidth=3)
        #    ax.plot(np.max(self.R_w)*np.cos(theta) , np.max(self.R_w)*np.sin(theta), 'k', linewidth=3)
        plt.show()
        
        
    def plot_Epitch(self):
        """
        Method to plot E pitch of ionisation
        """

        x=self.data_i['E']
        y=self.data_i['pitch']
            
        f=plt.figure()
        ax=f.add_subplot(111)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap, normed=True)
        f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('E')
        ax.set_ylabel('pitch')
        plt.show()
        
        
        
class fbm:
    """
    Manages the (tricky) output file for fbms from transp
    w3.pppl.gov/~pshare/transp/fbm.doc  
    """
    def __init__(self, infile_n):
        """
        R2D: X of MC zones
        Z2D: Z of MC zones
        plt.plot(R2D, Z2D) will give the picture of the 2D grid
        X2D: rho of MC zones
        TH2D: theta of MC zones
        A_D_NBI: pitch bin
        E_D_NBI: energy bins

        F_D_NBI: particle distribution function [#/cm^3/eV/d(omega/4pi)]        
        the actual distribution function f(E,p,R,Z) = 0.5*F_D_NBI (by integration
        over the spherical angle)
        
        """
        self.infile_n = infile_n
        self.infile = nc.Dataset(infile_n)
        self._read_info()        
        
        self.nbins=100        
        
        self.dict_dim = collections.OrderedDict([('R',[]),('z',[]),('pitch',[]),('E',[])])
        self.name_dim = {'R':'R2D', 'z':'Z2D', 'pitch':'A_D_NBI', 'E':'E_D_NBI'}
        self._read_dim()
        self.dict_dim['R'] *= 100.
        self.dict_dim['z'] *= 100.
        self.dict_dim_MCgrid = dict.copy(self.dict_dim)

        self.fdist_MCgrid = 0.5*self.infile.variables['F_D_NBI'][:]*1e6 #converting to m
        
        # This fdist_notnorm is the same as ascot distributions, this is useful because
        # I can copy the methods in this way       
        self.fdist_notnorm = np.zeros((self.shape_dim['E'], self.shape_dim['pitch'], \
                    self.nbins, self.nbins), dtype=float)        
        self._RZgrid()
        
        self.norm=1.
        
    def _read_info(self):
        """
        """         
        self.runid = ''
        self.tok = ''        
        runid = self.infile.variables['TRANSP_RUNID'][:]
        for i in runid:
            self.runid += i
        tok   = self.infile.variables['TOKAMAK'][:]
        for i in tok:
            self.tok += i
        self.time  = round(self.infile.variables['TIME'][:], 3)
        
    def _read_dim(self):
        """
        Hidden method to read the abscissae
        """
        self.shape_dim = self.dict_dim.copy()
        for i, dim in enumerate(self.dict_dim.iterkeys()):
            self.dict_dim[dim] = self.infile.variables[self.name_dim[dim]][:]
            self.shape_dim[dim] = np.size(self.dict_dim[dim])

    def _integrate_spacex(self, x, ax):
        """
        Hidden method to integrate over space and something else (pitch, E, mu...)
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        return np.trapz(self.f_space_int, self.dict_dim[x], axis=ax)    
            
    def _integrate_space(self):
        """
        Function to integrate over (R,z)
        """
        dist_toint = self.fdist_notnorm[:,:,:,:]/self.norm

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #pitch,E            


    def _RZgrid(self):
        """
        Converts from MC grid to R,Z common grid
        """
        print "Converting from MC grid to (R,Z)"
        x,y=self.dict_dim_MCgrid['R'], self.dict_dim_MCgrid['z']
        
        self.dict_dim['R'] = np.linspace(min(x), max(x), num=self.nbins, dtype=float)
        self.dict_dim['z'] = np.linspace(min(y), max(y), num=self.nbins, dtype=float)
        
        for i in range(self.shape_dim['pitch']):
            for j in range(self.shape_dim['E']):
                z = self._make_2d_fdist(self.dict_dim_MCgrid['R'], \
                                self.dict_dim_MCgrid['z'], self.fdist_MCgrid[:,i,j]) 
                self.fdist_notnorm[j,i,:,:] = z
             
    def _integrate_spacep(self):
        """
        hidden method to integrate over (space,pitch)
        """
        self.f_spacep_int = self._integrate_spacex('pitch', ax=1)        

    def _integrate_spaceE(self):
        """
        hidden method to integrate over (space,pitch)
        """
        self.f_spaceE_int = self._integrate_spacex('E', ax=0)    

    def plot_spaceE(self, norm):
        """
        plot 1D (pitch, int_space (int_E (fdist)))
        """
        try:
            self.f_spaceE_int.mean()
        except:
            self._integrate_spaceE()
        
        self.xplot = self.dict_dim['pitch']
        self.yplot = self.f_spaceE_int

        self._plot_1d('pitch', norm=norm)
        
    def plot_spacep(self, norm):
        """
        plot 1D (energy, int_space (int_pitch (fdist)))
        """
        try:
            self.f_spacep_int.mean()
        except:
            self._integrate_spacep()
        
        self.xplot = self.dict_dim['E']
        self.yplot = self.f_spacep_int
            
        self._plot_1d('Energy', norm=norm)


    def plot_Epitch(self):
        """
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        x,y = self.dict_dim['pitch'], self.dict_dim['E']
        title = self.runid + ' ' + str(self.time)
        self._plot_2d(x,y, self.f_space_int, r'$\xi$', r'E [eV]', title)

    def _make_2d_fdist(self, x,y,z):
        """
        Converts from MC irregular grid to 2D grid of distribution function
        """
        nbin_compl=complex(0,self.nbins)
        grid_x, grid_y = np.mgrid[min(x):max(x):nbin_compl, min(y):max(y):nbin_compl]
        grid_z = interp.griddata((x,y), z, (grid_x, grid_y), method='cubic', fill_value=0.)        
        return grid_z


    def _plot_1d(self, xlab, norm):
        """
        Hidden method to plot 1D functions
        """

        fig = plt.figure()
        ax  = fig.add_subplot(111)
        title = 'Distribution '
        if norm==1:
            self.yplot = self.yplot/self.norm
            title += ' NORM'
        ax.plot(self.xplot, self.yplot, linewidth=2.3)

        ax.set_xlabel(xlab)
        ax.set_title(title)
        
        plt.show()
        
    def _plot_2d(self,x,y,z, xlabel, ylabel, title):
        """
        """        
        fig = plt.figure()
        ax = fig.gca()
        CS = ax.contour(x,y,z, 20, cmap=my_cmap)
        plt.colorbar(CS)
        ax.set_xlabel(xlabel), ax.set_ylabel(ylabel)
        ax.set_title(title)
        plt.show()
        
class fbm_time:
    """
    Class as fbm including time dependence
    """
    def __init__(self, runid):
        self.runid = runid
        file_list = glob.glob('*fi*.cdf') # Get all the cdf in the current directory
        self.nslices=len(file_list)        
        self.timeslices = np.array([fbm for _ in range(self.nslices)])
        print "Number of timeslices: ", str(self.nslices)
        
        
        for i in range(self.nslices):
            fname = runid+'_fi_'+str(i+1)+'.cdf'
            self.nslices += 1
            print fname
            try:
                tmp = fbm(fname)
                self.timeslices[i] = tmp
            except:
                continue
        self._integratespace()
        
    def _integratespace(self):
        """
        """
        for i in self.timeslices:
            i._integrate_space()
            
    def plot_frames(self):
        """
        """
        
        n_cols = 3
        n_rows = self.nslices/3
        if self.nslices%3!=0:
            n_rows=n_rows+1
        vmax = np.max([i.f_space_int for i in self.timeslices])
        f, axarr = plt.subplots(n_cols, n_rows, sharex=True, sharey=True)        
        for i in range(self.nslices):
            ind_row = i/n_rows
            ind_col = i%n_rows
            CB=axarr[ind_row, ind_col].contourf(self.timeslices[i].dict_dim['pitch'], \
                    self.timeslices[i].dict_dim['E'], self.timeslices[i].f_space_int, \
                    20, cmap=my_cmap, vmin=0, vmax=vmax)
            #plt.colorbar(CB)        
            axarr[ind_row, ind_col].set_xlabel(r'$\xi$')
            axarr[ind_row, ind_col].set_ylabel(r'E [keV]')
            axarr[ind_row, ind_col].set_xlim([-1., 1.])
            axarr[ind_row, ind_col].set_ylim([0, 30000])
            axarr[ind_row, ind_col].set_title(str(self.timeslices[i].time))
        
        f.tight_layout()
        f.subplots_adjust(right=0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.colorbar(CB, cax=cbar_ax)        
        
        plt.show()
       
    def make_gif_Ep(self):
        """
        http://superfluoussextant.com/making-gifs-with-python.html
        """
        cwd = os.getcwd()
        tmpdir = cwd+'tmp_gifcreator/'
        self._make_Ep_singlefigures(tmpdir)
        os.chdir(tmpdir)
        
        gif_name = cwd+'Ep'+self.runid
        file_list = glob.glob('*.png') # Get all the pngs in the current directory
        #list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0])) # Sort the images by #, this may need to be tweaked for your use case

        with open('image_list.txt', 'w') as file:
            for item in file_list:
                file.write("%s\n" % item)
            
        os.system('convert -delay 100 @image_list.txt {}.gif'.format(gif_name)) # On windows convert is 'magick'
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        
    def _make_Ep_singlefigures(self, tmpdir):
        """
        """
        olddir = os.getcwd()
        os.makedirs(tmpdir); os.chdir(tmpdir)

        vmax = np.max([i.f_space_int for i in self.timeslices])
        for i, el in enumerate(self.timeslices):
            x,y = el.dict_dim['pitch'], el.dict_dim['E']
            z = el.f_space_int
            xlab, ylab = r'$\xi$', r'E [keV]'
            title = str(el.time)
            xlim, ylim = [-1., 1.], [0, 30000]
            self._save_singleframe(x, y, z, xlab, ylab, title, xlim, ylim, vmax, str(i))
        os.chdir(olddir)
        
        
    def _save_singleframe(self, x, y, z, xlab, ylab, title, xlim, ylim, vmax, savename):
        """
        """
        interactive(False) 
        levels = np.linspace(0,vmax,20)
        f = plt.figure()
        ax = f.add_subplot(111)        
        CB=ax.contourf(x,y,z, 20, cmap=my_cmap, vmin=0, vmax=vmax, levels=levels)
        plt.colorbar(CB)   
        ax.grid('on')
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title(title)
        plt.savefig(savename+'.png', bbox_inches='tight')
        interactive(True)