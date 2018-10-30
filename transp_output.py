#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, interactive, ticker
import math, collections
import scipy.interpolate as interp
import glob, os, shutil
from ascot_utils import common_style, limit_labels, _cumulative_plot, _plot_2d, _plot_1d

col = ['k','r','b','m','g','c']
col2=np.copy(col)
col2[0]='y'

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

class output_1d:
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
        self.runid = self.fname[-14:-4]
        self.run = self.runid[0:5]
        self._set_vars()
        self._calculate_all()
        
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
        self._kinetic_vars()
        self._nb_vars() 
        self._performance_vars()

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
        self._calculate_scalars()
        self._calculate_eff()
        self._calculate_Ec()
        self._calculate_gi()
        self._calculate_tau()
        
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
        self.darea = np.array(self.file.variables['DAREA'][:])*1e-4
        self.dvol  = np.array(self.file.variables['DVOL'][:])*1e-6
        self.rho = np.linspace(0, 1, num=len(self.darea[0,:]), dtype=float)

        self._imp_vars()
        
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
        self.imp_names, self.imp_vars = \
            self._fill_dict(keys, varnames)        
        
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
            self._fill_dict(keys, varnames)
        self.kin_vars['ne'] *= 1e6
        self.kin_vars['ni'] *= 1e6

    def _nb_vars(self):
        """ NB-related variables

        Gathers: 
        - INPUT: beam energy, beam power, beam species, beam geometry
        - simulation results (nfast, orbloss, cx, pi, pe, pth, nbcd, unbcd, 
                              ptot(ini-end), press of EP) from netcdf file
        Parameters:
            None
        Attributes:
            None
        Note:
            Everything is imperial system (i.e. curb is A/cm**2, n is cm**-3, etc.)

            If two D beams are present, usually the first is the heating beam

        """
        keys     = ['P', 'E']
        varnames = ['PBINJ_D', 'EINJAV_D' ]
        self.nb_in_names = dict.fromkeys(keys)
        self.nb_in_vars  = dict.fromkeys(keys)
        self.nb_in_names, self.nb_in_vars = \
            self._fill_dict(keys, varnames)        
        self.inj_index = np.where(self.nb_in_vars['P']>0.9*np.max(self.nb_in_vars['P']>0.9))[0]

        keys     = ['st']
        varnames = ['PBSHINE_D']
        # self.nb_ioniz_names = dict.fromkeys(keys)
        # self.nb_ioniz_vars  = dict.fromkeys(keys)
        self.nb_ioniz_names, self.nb_ioniz_vars = \
            self._fill_dict(keys, varnames)         
        
        keys     = ['n','orbloss', 'cx1', 'cx2', 'pi', 'pe',\
                    'th', 'nbcd', 'unbcd', 'pMC', 'press_EP']
        keys_red = ['pi', 'pe', 'th']
        varnames = ['BDENS_D','BPLIM_D', 'BPCXX_D', 'BPCXI_D', 'PBI_D', 'PBE_D', \
                    'PBTH_D', 'CURB', 'UCURB', 'PBDEPMC_D', 'PTOWB']
        varnames1 = ['PBI01_TOT', 'PBE01_TOT', 'PBTH01']
        varnames2 = ['PBI02_TOT', 'PBE02_TOT', 'PBTH02']
        # self.nb_FKP_names = dict.fromkeys(keys)
        # self.nb_FKP_vars  = dict.fromkeys(keys)
        # self.Hnb_FKP_names = dict.fromkeys(keys_red)
        # self.Hnb_FKP_vars  = dict.fromkeys(keys_red)
        # self.Dnb_FKP_names = dict.fromkeys(keys_red)
        # self.Dnb_FKP_vars  = dict.fromkeys(keys_red)

        self.nb_FKP_names, self.nb_FKP_vars = self._fill_dict(keys, varnames)
        # check if needed to split in the two beam components
        if varnames2[0] in self.file.variables.keys():
            self.Hnb_FKP_names, self.Hnb_FKP_vars = self._fill_dict(keys_red, varnames1)
            self.Dnb_FKP_names, self.Dnb_FKP_vars = self._fill_dict(keys_red, varnames2) 
            #self.Hnb_FKP_vars['n'] *= 1e6
            #self.Dnb_FKP_vars['n'] *= 1e6
            self.Hnb_FKP_vars['pi'] = self.Hnb_FKP_vars['pi']*1.e6
            self.Hnb_FKP_vars['pe'] = self.Hnb_FKP_vars['pe']*1.e6
            self.Dnb_FKP_vars['pi'] = self.Dnb_FKP_vars['pi']*1.e6
            self.Dnb_FKP_vars['pe'] = self.Dnb_FKP_vars['pe']*1.e6
            
        # conversion to metric system
        self.nb_FKP_vars['n'] *= 1e6

        self.nb_FKP_vars['nbcd']*=1e4
        self.nb_FKP_vars['unbcd']*=1e4

        self.nb_FKP_vars['pi'] *= 1e6
        self.nb_FKP_vars['pe'] *= 1e6

       
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
        keys     = ['jboot', 'jbootsau', 'jbootneogk', 'jbootneo', 'jbootsauor']
        varnames = ['CURBS', 'CURBSSAU', 'CURBSNEO', 'CURBSWNC', 'CURBSSAU0']
        self.perf_names = dict.fromkeys(keys)
        self.perf_vars  = dict.fromkeys(keys)
        self.perf_names, self.perf_vars = \
            self._fill_dict(keys, varnames)
        for el in self.perf_vars:
            self.perf_vars[el] *= 1e4


    def _fill_dict(self, keys, varnames):
        """Fills dictionaries

        Given the name of the keys and the name of the variables,
        fills the dictionaries

        Parameters:
            keys      (array) : keys of the output dictionaries
            varnames  (array) : name of NETCDF file variables
        Attributes:
            name_dict (dict) : dict associating keys with varnames
            var_dict  (dict) : dict associating actual variables
        Note:

        """
        name_dict = dict.fromkeys(keys)
        var_dict  = dict.fromkeys(keys)
        
        for ii,kk in enumerate(keys):
            name_dict[kk] = varnames[ii]
            tmpdata = self.file.variables[name_dict[kk]][:]
            var_dict[kk] = tmpdata
            
        return name_dict, var_dict
    
    
    def _calculate_Ec(self):
        """ Critical energy

        Computes critical energy in the plasma
        Ec = 14.8*Te*(A^{3/2}/ne*sum(nj*Zj/Aj))^2/3

        Parameters:
            None
        Attributes:
            None    
        Note:

        """
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
            
    def _calculate_tau(self):
        """ tau_s

        Computes EP slowing-down time
        tau_s=ts/3*ln(1+(E/EC)**1.5)
        ts=6.27e8*A Te[eV]**1.5/(Z**2 ne(cm^-3) lnLambda)
        valid for Z=1

        Parameters:
            None
        Attributes:
            None    
        Note:

        """
        try:
            self.Ec.mean()
        except:
            self._calculate_Ec()

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
            taus=ts/3.*np.log((1+(E0/Ec[it,:])**1.5))
            self.taus[it,:] = taus
            self.taus_mean[it] = np.trapz(taus, self.rho)

        
    def _calculate_gi(self):
        """ Gi

        Computes EP deposition to ions
        Gi = Ec/E0 int[0, E0/Ec]dy/(1+y**1.5)

        Parameters:
            None
        Attributes:
            None    
        Note:

        """

        try:
            self.Ec.mean()
        except:
            self._calculate_Ec()
        Ec = self.Ec; E0 = self.nb_in_vars['E']
        self.gi = np.zeros(np.shape(Ec), dtype=float)
        self.gi_vavg = np.zeros(self.nt, dtype=float)
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
            self.gi_vavg[time_ind] = np.trapz(self.gi[time_ind,:], self.rho)   
        self.gi_mean = np.mean(self.gi_vavg[self.inj_index])

    def _time_to_injind(self, time=[0]):
        """ times to indexes

        Gets the indexes of injection time relative to the times requested

        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            ind (arr)  : array with the indexes where to plot the lines
        Note:

        """
        if time[0]!=0:
            ind=[]
            for t in time:
                ind.append(np.argmin(self.t[self.inj_index]-t<0.))
            ind=np.array(ind,dtype=int)
        else:
            if self.nt > 5:
                ind = np.array([0])
                #ind = np.linspace(np.min(self.inj_index), np.max(self.inj_index)-1, 5)
                ind = ind.astype(int)
        return ind

    def _time_to_ind(self, time=[0]):
        """ times to indexes

        Gets the indexes relative to the times requested

        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            ind (arr)  : array with the indexes where to plot the lines
        Note:

        """
        if time[0]!=0:
            ind=[]
            for t in time:
                ind.append(np.argmin(self.t-t<0.))
            ind=np.array(ind,dtype=int)
        else:
            if self.nt > 5:
                ind = np.array([0])
                #ind = np.linspace(0, len(self.t)-1, 5)
                ind = ind.astype(int)

        return ind

    def plot_input(self, time=[0]):
        """ plot all the input

        Plots all the input parameters (vol-avgd and profiles)
        calling the two relative functions

        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            None    
        Note:
        """
        ind = self._time_to_ind(time)
        self.plot_input_1d(ind)
        self.plot_input_prof(ind)
        
        
    def plot_input_1d(self, time=[0], axne=0, axTe=0, ls='-'):
        """ plot 1d input

        Plots ne,ni,te,ti vol-avgd

        Parameters:
            time (arr) : array with the times where to plot the lines
        Attributes:
            None    
        Note:
        """
        common_style()

        if axne==0:
            f = plt.figure(figsize=(8, 8))
            axne = f.add_subplot(211)
            axTe = f.add_subplot(212, sharex=axne)
            fig_flag = 0 # this means it wasn't initialized 
        #axzf = f.add_subplot(313, sharex=axne)
        else:
            f=plt.gcf()
            fig_flag=1

        t=self.t
        for y, c, l in zip([self.ne_vavg, self.ni_vavg], ['k','r'], [r'e', r'i']):
            _plot_1d(t, y, ax=axne, color=c, label=l, ls=ls)
        for y, c, l in zip([self.Te_vavg, self.Ti_vavg], ['k','r'], [r'e', r'i']):
            _plot_1d(t, y, ax=axTe, color=c, label=l, ls=ls)

        if fig_flag==0:
            #====================================================
            # Correct ticks and xylabels
            #====================================================
            limit_labels(axne, r'Time [s]', r'$\langle n \rangle$ [$1/m^3$]','' )
            limit_labels(axTe, r'Time [s]', r'$\langle T \rangle$ [$eV$]','' )

        #====================================================
        # Plot vertical lines
        #====================================================
        ind = self._time_to_ind(time)
        for ax in [axne, axTe]:#, axzf]:
            if len(ind)!=1:
                for i, el in enumerate(ind):
                    ax.axvline(x=self.t[el], color=col[i], lw=2., linestyle='--')
            if fig_flag==0:
                ax.legend(loc='best')

        f.tight_layout()
        plt.show()    

       
    def plot_input_prof(self, time=[0]):
        """ plot 1d input

        Plots input quantities, such as ne, ni, Te, Ti profiles
        
        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            None    
        Note:
        """        
        common_style()
       
        f = plt.figure(figsize=(12,10))
        axne = f.add_subplot(221)
        axni = f.add_subplot(222, sharey=axne)
        axTe = f.add_subplot(223)
        axTi = f.add_subplot(224, sharey=axTe)
        ind = self._time_to_ind(time)

        for i, el in enumerate(ind):
            lab = r't = {:.2f}'.format(self.t[el])
            axne.plot(self.rho, self.kin_vars['ne'][el,:], col[i], label=lab, lw=2.)
            axni.plot(self.rho, self.kin_vars['ni'][el,:], col[i], label=lab, lw=2.)
            axTe.plot(self.rho, self.kin_vars['te'][el,:]*1e-3, col[i], label=lab, lw=2.)
            axTi.plot(self.rho, self.kin_vars['ti'][el,:]*1e-3, col[i], label=lab, lw =2.)

        #========================================================
        # SET TICK LOCATION
        #========================================================
        limit_labels(axne, r'$\rho$', r'$n_e$ [$1/m^3$]','' )
        limit_labels(axTe, r'$\rho$', r'$T_e$ [$keV$]','' )
        limit_labels(axni, r'$\rho$', r'$n_i$ [$1/m^3$]','' )
        limit_labels(axTi, r'$\rho$', r'$T_i$ [$keV$]','' )


        for ax in [axne, axni, axTe, axTi]:
            # Create your ticker object with M ticks
            M = 4
            yticks = ticker.MaxNLocator(M)
            xticks = ticker.MaxNLocator(M)
            # Set the yaxis major locator using your ticker object. You can also choose the minor
            # tick positions with set_minor_locator.
            ax.yaxis.set_major_locator(yticks)
            #ax.yaxis.set_minor_locator(yticks_m)
            ax.xaxis.set_major_locator(xticks)
            #========================================================

        #axne.legend(bbox_to_anchor=(0, .1, -0.1, .102))#, axTe.legend(loc='best')
        axTi.legend(loc='best')
        #f.text(0.01, 0.01, self.fname)

        f.tight_layout()
        #f.subplots_adjust(left=0.2)
        plt.show()    

    def plot_deposition(self, time=[0]):
        """
        """

        if len(time)==1:
            if self.nt > 5:
                ind = np.linspace(np.min(self.inj_index), np.max(self.inj_index)-1, 5)
                ind = ind.astype(int)
        else:
            ind=[]
            for t in time:
                ind.append(np.argmin(self.t[self.inj_index]-t<0.))
                print t, self.t[self.inj_index][ind]
            ind=np.array(ind,dtype=int)

        self._plot_deposition_1d(ind)
        self._plot_deposition_prof(ind)

    def plot_deposition_timeslice(self, timeslice):
        """
        """
        ind=np.argmin(self.t-timeslice<0)
        ind=np.array(ind)
        self._plot_deposition_prof(ind=[ind])
        

    def _plot_deposition_1d(self, time=[0], axp=0, axf=0, ls='-'):
        """
        """
        try:
            self.gi.mean()
        except:
            self._calculate_gi()
        try:
            self.pe.mean()
        except:
            self._calculate_scalars()
        
        common_style()
        ind = self._time_to_injind(time)
        if axp==0:
            f = plt.figure(figsize=(10, 8))
            axp = f.add_subplot(211)
            axf = f.add_subplot(212, sharex=axp)
            fig_flag = 0 # this means it wasn't initialized 
        else:
            f=plt.gcf()
            fig_flag=1            
        t=self.t; ptot=self.pi+self.pe
        for y, c, l in zip([self.pe, self.pi], ['k','r'], [r'e', r'i']):
            _plot_1d(t, y, ax=axp, color=c, label=l, ls=ls)
        for y, c, l in zip([self.pe, self.pi], ['k','r'], [r'e', r'i']):
            _plot_1d(t, y*100./ptot, ax=axf, color=c, label=l, ls=ls)

        axf.plot(self.t, self.gi_vavg*1e2, 'c', lw=2.3, label=r'G')

        if fig_flag==0:
            #====================================================
            # Correct ticks and xylabels
            #====================================================
            limit_labels(axp,r'Time [s]',r'P[kW]','')
            limit_labels(axf,r'Time [s]',r'%','')

        #====================================================
        # Plot vertical lines
        #====================================================
        for ax in [axp, axf]:
            if ind[0]!=0:
                for i, el in enumerate(ind):
                    ax.axvline(x=self.t[self.inj_index[el]], color=col[i], lw=2., linestyle='--')

        axp.set_ylim([0,160]); axf.set_ylim([0,100.])
        f.tight_layout()
        plt.show()     
        
        
    def _plot_deposition_prof(self, ind=[0], **kwargs):
        """
        Plots deposition to ions, electrons
        """
        common_style()
        try:
            self.nb_FKP_vars['pi'].mean()
        except:
            print "No Pi data"
            return
        if len(ind)==1:
            ind = np.linspace(0, len(self.inj_index)-1, 5, dtype=int)
        f = plt.figure(figsize=(18,6))
        axe = f.add_subplot(131)
        axi = f.add_subplot(132, sharex=axe, sharey=axe)
        axn = f.add_subplot(133, sharex=axe)
        #axj = f.add_subplot(144, sharex=axe)
        x=self.rho
        for index, i in enumerate(self.inj_index[ind]):
            lab=r't = {:.2f} s'.format(self.t[i])
            axi.plot(x, self.nb_FKP_vars['pi'][i,:]*1e-3, col[index], lw=2, label=lab)
            axe.plot(x, self.nb_FKP_vars['pe'][i,:]*1e-3, col[index], lw=2, label=lab)
            axn.plot(x, self.nb_FKP_vars['n'][i, :]/self.kin_vars['ne'][i,:]*100., col[index], lw=2, label=lab)
            #axj.plot(x, self.nb_FKP_vars['nbcd'][i,:]*1e-3, col[index], lw=2, label=lab)
        limit_labels(axe,r'$\rho$', r'$P_e$ [$kW/m^3$]','')
        limit_labels(axi,r'$\rho$', r'$P_i$ [$kW/m^3$]','')
        limit_labels(axn,r'$\rho$', r'$n_f/n_e$ [%]','')
        #limit_labels(axj,r'$\rho$', r'j shielded [$kA/m^2$]','')
        
        axn.legend(bbox_to_anchor=(1.05, 0.65), loc=2)
        #f.text(0.01, 0.01, self.fname)

        f.tight_layout()
        f.subplots_adjust(right=0.8)
        plt.show()   
        
        
    def plot_NBCD(self, time=[0]):
        """Plot current drive
        
        Makes a plot with NBCD (shielded and not), the shielding factor and the efficiency

        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            None    
        Note:        
        """
        common_style()

        plt.rc('legend', fontsize=15)

        try:
            self.nb_FKP_vars['nbcd'].mean()
        except:
            print "No nbcd data"
            return

        ind = self._time_to_injind(time)
        f = plt.figure(figsize=(10,8))
        axj  = f.add_subplot(224)
        axcd = f.add_subplot(221)
        axsh = f.add_subplot(222)
        axeff = f.add_subplot(223)

        for ax in [axcd, axsh, axeff]:
            if time[0]!=0:
                for i, el in enumerate(ind):
                    ax.axvline(x=self.t[self.inj_index[el]], color=col[i], lw=1.8, linestyle='--')

        axcd.plot(self.t[self.inj_index], self.nbcd[self.inj_index]*1e-3, 'k', lw=2.3, label='SH.')
        axcd.plot(self.t[self.inj_index], self.unbcd[self.inj_index]*1e-3, 'k:', lw=2.3, label='UNSH.')   
        axsh.plot(self.t[self.inj_index], 1.-self.shield[self.inj_index], 'k', lw=2.3)
        axeff.plot(self.t[self.inj_index], self.eff[self.inj_index]*100., 'k', lw=2.3)
        for i, el in enumerate(ind):
            axj.plot(self.rho, self.nb_FKP_vars['nbcd'][el,:]*1e-3, \
                     col[i],lw=2.3,label=r't = {:.2f} s'.format(self.t[el]))
        limit_labels(axj, r'$\rho$', r'$j^{SH}$ [$kA/m^2$]')
        limit_labels(axcd, r't [s]', r'$I_{CD}$ [$kA$]')
        axcd.legend(loc='best')
        limit_labels(axsh, r't [s]', r'$1-I_{SH}/I_{UN}$')
        limit_labels(axeff, r't [s]', r' $\eta \left[\frac{10^{18} A}{W m^2}\right]$', r'NBCD efficiency')
        axeff.set_ylim([0,1.])
        axj.legend(bbox_to_anchor=(-2.35,1.65), loc=2)
        #f.text(0.01, 0.01, self.runid)
        f.tight_layout();
        f.subplots_adjust(left=0.22, right=0.9)

        plt.show()

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
        self.nf_vavg = np.zeros(self.nt, dtype=float)
        
        for i in range(self.nt):
            dvol = self.dvol[i,:]
            self.ne_vavg[i] = np.dot(self.kin_vars['ne'][i,:],dvol)/np.sum(dvol)
            self.ni_vavg[i] = np.dot(self.kin_vars['ni'][i,:],dvol)/np.sum(dvol)
            self.Te_vavg[i] = np.dot(self.kin_vars['te'][i,:],dvol)/np.sum(dvol)
            self.Ti_vavg[i] = np.dot(self.kin_vars['ti'][i,:],dvol)/np.sum(dvol)
            self.nf_vavg[i] = np.dot(self.nb_FKP_vars['n'][i,:],dvol)/np.sum(dvol)
            
        self.ne_mean = np.mean(self.ne_vavg[self.inj_index])
        self.ni_mean = np.mean(self.ni_vavg[self.inj_index])
        self.Te_mean = np.mean(self.Te_vavg[self.inj_index])
        self.Ti_mean = np.mean(self.Ti_vavg[self.inj_index])
        self.nf_mean = np.mean(self.nf_vavg[self.inj_index])
        
    def _calculate_scalars(self):
        """ Scalar calculations

        Calculates scalar variables (integrals of profiles)

        Parameters:
            None
        Attributes:
            None    
        Note:
        """

        self.pe = np.zeros(len(self.t), dtype=float)
        self.pi = np.zeros(len(self.t), dtype=float)
        self.pth = np.zeros(len(self.t), dtype=float)
        #difference between H and D beam
        self.Hpe = np.zeros(len(self.t), dtype=float)
        self.Hpi = np.zeros(len(self.t), dtype=float)
        self.Hpth = np.zeros(len(self.t), dtype=float)
        self.Dpe = np.zeros(len(self.t), dtype=float)
        self.Dpi = np.zeros(len(self.t), dtype=float)
        self.Dpth = np.zeros(len(self.t), dtype=float)

        self.nbcd = np.zeros(len(self.t), dtype=float)
        self.unbcd = np.zeros(len(self.t), dtype=float)
        self.jboot = np.zeros(len(self.t), dtype=float)
        self.shield = np.zeros(len(self.t), dtype=float)        

        self.jbootSAU   = np.zeros(len(self.t), dtype=float)
        self.jbootneogk = np.zeros(len(self.t), dtype=float)
        self.jbootneo   = np.zeros(len(self.t), dtype=float)
        self.jbootSAUor = np.zeros(len(self.t), dtype=float)

        for i in range(len(self.t)):
            self.pe[i] = np.dot(self.nb_FKP_vars['pe'][i,:], self.dvol[i,:])
            self.pi[i] = np.dot(self.nb_FKP_vars['pi'][i,:], self.dvol[i,:])
            self.pth[i] = np.dot(self.nb_FKP_vars['th'][i,:], self.dvol[i,:])
            if 'PBI02_TOT' in self.file.variables.keys():
                self.Hpe[i] = np.dot(self.Hnb_FKP_vars['pe'][i,:], self.dvol[i,:])
                self.Hpi[i] = np.dot(self.Hnb_FKP_vars['pi'][i,:], self.dvol[i,:])
                self.Hpth[i] = np.dot(self.Hnb_FKP_vars['th'][i,:], self.dvol[i,:])
                self.Dpe[i] = np.dot(self.Dnb_FKP_vars['pe'][i,:], self.dvol[i,:])
                self.Dpi[i] = np.dot(self.Dnb_FKP_vars['pi'][i,:], self.dvol[i,:])
                self.Dpth[i] = np.dot(self.Dnb_FKP_vars['th'][i,:], self.dvol[i,:])

            self.nbcd[i] = np.dot(self.nb_FKP_vars['nbcd'][i,:], self.darea[i,:])
            self.unbcd[i] = np.dot(self.nb_FKP_vars['unbcd'][i,:], self.darea[i,:])            
            self.jboot[i] = np.dot(self.perf_vars['jboot'][i,:], self.darea[i,:])
            self.shield[i] = self.nbcd[i]/self.unbcd[i]           

            self.jbootSAU[i]   = np.dot(self.perf_vars['jbootsau'][i,:], self.darea[i,:])
            self.jbootneogk[i] = np.dot(self.perf_vars['jbootneogk'][i,:], self.darea[i,:])
            self.jbootneo[i]   = np.dot(self.perf_vars['jbootneo'][i,:], self.darea[i,:])
            self.jbootSAUor[i] = np.dot(self.perf_vars['jbootsauor'][i,:], self.darea[i,:])


        self.pcx = self.nb_FKP_vars['cx1']+self.nb_FKP_vars['cx2']
        self.pol = self.nb_FKP_vars['orbloss']
        self.psh = self.nb_ioniz_vars['st']

        ind = self.inj_index
        self.psh_mean = np.mean(self.psh[ind])
        self.pcx_mean = np.mean(self.pcx[ind])
        self.pol_mean = np.mean(self.pol[ind])
        self.pin_mean = np.mean(self.nb_in_vars['P'][ind])
        self.pi_mean  = np.mean(self.pi[ind])
        self.pe_mean  = np.mean(self.pe[ind])
        self.pth_mean = np.mean(self.pth[ind])
        self.pabs_mean= self.pi_mean+self.pe_mean+self.pth_mean 
        self.nbcd_mean = np.mean(self.nbcd[ind])        
        
    def print_power_balance(self):
        """ power balance

        Prints the average power balance

        Parameters:
            None
        Attributes:
            None    
        Note:
        """            
        pbal = np.mean(self.nb_in_vars['P']-self.psh - \
            self.pe-self.pi-self.pth-self.pol-self.pcx)*1e-6
        print "POWER BALANCE:", pbal
        print "Pbalance/Pin ", pbal/np.mean(self.nb_in_vars['P'])*1e6
        print "psh+pol ", np.mean(self.pol+self.psh)*1e-6
        print "pe+pi+pth ", np.mean(self.pe+self.pi+self.pth)*1e-6
        print "pcx " , np.mean(self.pcxprof)*1e-6

        
    def plot_performance(self):
        """
        """
        common_style()
        try:
            self.perf_vars['jboot'].mean()
        except:
            print "No bootstrap data"
            return
        self._calculate_scalars()

        f = plt.figure()
        axcd = f.add_subplot(111)          
        axcd.plot(self.t, self.jboot*1e-3, 'k', lw=2.3, label='JBS')
        axcd.plot(self.t, self.jbootSAU*1e-3, 'r', lw=2.3, label='SAUTER')
        axcd.plot(self.t, self.jbootneogk*1e-3, 'b', lw=2.3, label='NEO GK')
        axcd.plot(self.t, self.jbootneo*1e-3, 'c', lw=2.3, label='NEO')
        axcd.plot(self.t, self.jbootSAUor*1e-3, 'm', lw=2.3, label='SAUTER OR')


        limit_labels(axcd, r't [s]', r'Bootstrap current [kA]')
        f.tight_layout(), f.suptitle(self.fname)
        axcd.legend(loc='best')
        axcd.grid('on')
        plt.show()


    def plot_powerbalance(self, time=[0], ls='-'):
        """ plot pow balance

        Plots the power balance with cumulative plot
       
        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            None    
        Note:
        """   
        #################################
        # defining vars to plot
        #################################
        ind = self._time_to_ind(time)
        x,y,yH,yfrac = self._init_pbalance_vars()

        labels = ['Abs.', 'CX', 'O. L.', 'Shine-thr.']
        xlabel = r'Time [s]'
        ylabel = r'Power [MW]'
        
        #################################
        # Plot
        #################################
        f=plt.figure(figsize=(15,5))
        axfrac = f.add_subplot(122)
        self._plot_pbalance_frac(x, yfrac, axfrac, ls=ls)
        axabs = f.add_subplot(121)
        _cumulative_plot(x,y,labels, xlabel, ylabel, col2, ax=axabs, title='')

        if len(time)!=1:
            for i, el in enumerate(ind):
                axfrac.axvline(x=self.t[self.inj_index[el]], color='k', lw=2.3, linestyle='--')       
                axabs.axvline(x=self.t[self.inj_index[el]], color='k', lw=2.3, linestyle='--')      
        f.subplots_adjust(right=0.8)        
        f.tight_layout()

    def _init_pbalance_vars(self):
        """initialize vars for pbalance
        
        Parameters:
        Attributes:
            None    
        Note:        
        
        """
        x = self.t
        y = [(self.pi+self.pe+self.pth)*1e-6, self.pcx*1e-6, self.pol*1e-6, self.psh*1e-6]
        yH = [(self.Hpi+self.Hpe+self.Hpth)*1e-6, self.pcx*1e-6, self.pol*1e-6, self.psh*1e-6]
        pin = self.nb_in_vars['P']
        yfrac = y/pin*1e6
        return x,y,yH,yfrac

    def _plot_pbalance_frac(self, x, y, ax,  ylim=[0,50], ls='-', leg=1):
        """ plot frac pow balance

        Plots the fractional power balance with cumulative plot
       
        Parameters:
            x (arr): time values on abscissa
            y (arr): set of y values (n x times). will plot all elements in y
            ax (obj): axis object where to plot
            ylim (array): limits on y. default to [0,50]
            ls : style of the line. default to '-'
        Attributes:
            None    
        Note:
        """ 
        labels = ['Abs.', 'CX', 'O. L.', 'Shine-thr.']

        for i, el in enumerate(y):
            ax.plot(x, el*100., color=col2[i], lw=2.3, label=labels[i], linestyle=ls)
        ax.set_ylim(ylim)
        #f.text(0.01, 0.01, self.runid)
        if leg==1:
            limit_labels(ax, r'Time [s]', r'Fraction [%]')        
            #ax.legend(loc=2, bbox_to_anchor=(1.05, .7))
            ax.legend(loc='best')

    def _calculate_eff(self):
        """ calculate efficiency

        Method to calculate the current-drive efficiency from the beam:
        eta = R0*n_e*I_CD/P [10^20 A / (W m^2)]
       
        Parameters:
            None
        Attributes:
            None    
        Note:
        """                
        try:
            self.nbcd.mean()
        except:
            self._calculate_scalars()
        #Computing the power coupled to plasma: Pini-Pend+Pres
        P = self.nb_in_vars['P']
        #ind = P>0.
        R0 = self.file.variables['RAXIS'][:]*0.01
        # average density in 10^20
        ne_avg = self.ne_vavg[:]*1e-20
        #############################
        self.eff = R0*ne_avg*np.abs(self.nbcd)/P
        self.eff[np.isinf(self.eff)] = 0.
        self.eff[np.isnan(self.eff)] = 0.
        self.eff_mean = np.mean(self.eff[self.inj_index])
#        print "R0 ", R0, 'm'
#        print "ne avg", ne_avg, '10e20 m^-3'
#        print "Ip ", self.nbcd, 'A'
#        print "P  ", P, 'W'
#        print "CD efficiency: ", self.eff, " [10^20 A / (W m^2)]"

    def plot_radar_plost(self):
        """
        """
        self._calculate_scalars()
        N=3
        values = [self.psh_mean/self.pin_mean*100., self.pol_mean/self.pin_mean*100., \
                    self.pcx_mean/self.pin_mean*100.]

        categories = ['ST', 'OL', 'CX']
        _radar_plot(N, values, categories,  title='Losses')
        
        
    def plot_radar_deposition(self):
        """
        radar plot of absorbed power, NBCD, nf/ne, %powers to ions, efficiency, momentum        
        """
        N=4
        
        values = [self.pabs_mean/self.pin_mean*100., self.nbcd_mean*1e-3, \
                    self.nf_mean/self.ne_mean*100., self.gi_mean*100.]#, \
                    #self.eff_mean ]
        categories = [r'$P_{ABS}$', r'$I_{CD}$', r'$n_f/n_e$ [%]', r'$G_i$']#, r'$\eta$']
        _radar_plot(N, values, categories, title='Deposition')


    def write_NB_out(self, timeslice):
        """ write NB to ascii

        Method writing a timeslice of a simulation
        rho, power transferred to electrons and ions, the fast ion density and pressure and the current drive

        Parameters:
            None
        Attributes:
            None    
        Note:
        """
        ind = np.argmin(np.abs(self.t-timeslice))
        data = [self.rho, self.nb_FKP_vars['pe'][ind,:], self.nb_FKP_vars['pi'][ind,:],\
                self.nb_FKP_vars['n'][ind,:], self.nb_FKP_vars['press_EP'][ind,:],\
                self.nb_FKP_vars['nbcd'][ind,:]]
        header = 'rho tor\t  Pe[W/m3]\t Pi[W/m3]\t n[1/m3]\t press[Pascal]\t nbcd[A/m2]'
        out_fname='output_'+str(timeslice*100.)+'.dat'
        np.savetxt(out_fname, np.transpose(data), fmt='%.5e', delimiter='\t', header=header)

    def write_prof_in(self, timeslice):
        """ write profiles to ascii

        Method writing a timeslice of a simulation
        rho, ne, ni, Te, Ti, Zeff

        Parameters:
            None
        Attributes:
            None    
        Note:
        """
        ind = np.argmin(np.abs(self.t-timeslice))
        data = [self.rho, self.kin_vars['ne'][ind,:], self.kin_vars['ni'][ind,:], self.kin_vars['te'][ind,:], self.kin_vars['ti'][ind,:]]
        header = 'rho tor\t  ne[1/m3]\t ni[1/m3]\t Te[eV]\t\t Ti[eV]'
        out_fname='input_'+str(timeslice*100.)+'.dat'
        np.savetxt(out_fname, np.transpose(data), fmt='%.5e', delimiter='\t', header=header)

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
        self.n0_tot = tot_source

        common_style()
        f=plt.figure(); ax=f.add_subplot(111)
        ax.plot(self.t, tot_source[:, -1]*1e6, 'k', lw=2.3, label='tot')
        ax.plot(self.t, v_source[:,-1]*1e6, 'b', lw=2.3, label='Volume')
        ax.plot(self.t, halob[:,-1]*1e6, 'm', lw=2.3, label='halo')
        ax.plot(self.t, first_fastn[:,-1]*1e6, 'c', lw=2.3, label='fastn')
        ax.plot(self.t, w_source[:,-1]*1e6, 'r', lw=2.3, label='wall')
        ax.plot(self.t, Drecy[:,-1]*1e6+Dnrecy[:,-1]*1e6+Drec[:,-1]*1e6, 'g', lw=2.3, label='Recy')
        ax.plot(self.t, CXfastn[:,-1]*1e6,'y', lw=2.3, label='CX')
        ax.set_xlabel(r't [s]'); ax.set_ylabel(r'n0 [1/m3]')
        ax.legend(loc='best'); ax.grid('on')
        plt.show()

    def _n0edge(self):
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

        CXfastn = self.file.variables['N0BCXD0'][:]
        first_fastn = self.file.variables['N0BD0'][:]
        halob = self.file.variables['N0BH_D'][:]
        Drecy = self.file.variables['N0RC_D_D'][:]
        Dflow= self.file.variables['N0GF_D_D'][:]
        Dsflow= self.file.variables['N0SGF_D'][:]
        Dnrecy= self.file.variables['N0SRC_D'][:]
        Drec= self.file.variables['N0V0_D'][:]
         
        tot_source = v_source+w_source+first_fastn+CXfastn+Drecy+Dflow+halob+Dsflow+Dnrecy+Drec
        self.n0_tot = tot_source


    def _calculate_wmhd(self):
        """
        Calculate the perpendicular stored energy
        wi, we, wimp = perpendicular energies of bulk ions, electrons and impurities
        wth= only thermal component of the energy
        wprp = perpendicular component of fast ion energy
        """
        if self.flag_beam == 1:
            wi  = np.multiply(self.kin_vars['ni']-self.nb_FKP_vars['n'],self.kin_vars['ti'])*1.602e-19*1.5
        else:
            wi  = np.multiply(self.kin_vars['ni'],self.kin_vars['ti'])*1.602e-19*1.5
        we  = np.multiply(self.kin_vars['ne'],self.kin_vars['te'])*1.602e-19*1.5
        wim = np.multiply(self.imp_vars['nimp'],self.kin_vars['ti'])*1.602e-19*1.5
        wth = we+wi+wim
        self.wth=np.zeros(len(self.t)); self.wtot=np.copy(self.wth); self.wprp=np.copy(self.wth)
        
        for i in range(len(self.t)):
            self.wth[i]  = np.dot(wth[i], self.dvol[i,:])
            if self.flag_beam==1:
                self.wprp[i] = np.dot(self.file.variables['UBPRP'][i], self.dvol[i,:])
            else:
                self.wprp[i]=0.
        self.wtot=self.wth+self.wprp
