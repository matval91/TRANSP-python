#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
from __future__ import print_function

import numpy as np

from pytransp.classes.transp_output import transp_output
import pytransp.trutils.transp_utils as tu

import utils.plot_utils as au
import pytransp.plot.plot_nb as plot_nb

class transp_heating(transp_output):
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
        transp_output.__init__(self, fname)
        if 'BDENS' in self.file.variables.keys():
            self._nb_vars() 
            self._calculate_all_nb()
            self.flag_beam=1
    
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
        self.nb_in_names, self.nb_in_vars = \
            tu._fill_dict(self.file, keys, varnames)        
        self.inj_index = np.where(self.nb_in_vars['P']>0.9*np.max(self.nb_in_vars['P']>0.9))[0]

        keys     = ['st']
        varnames = ['PBSHINE_D']
        if 'PBSHINE_H' in self.file.variables.keys():
            keys.append('st_H')
            varnames.append('PBSHINE_H')
        self.nb_ioniz_names, self.nb_ioniz_vars = \
            tu._fill_dict(self.file, keys, varnames)   
            
        # GLOBAL FAST ION VARIABLES
        keys     = ['n','orbloss', 'cx1', 'cx2', 'pi', 'pe',\
                    'th', 'nbcd', 'unbcd', 'press_EP']
        varnames = ['BDENS','BPLIM', 'BPCXX', 'BPCXI', 'PBI', 'PBE', \
                    'PBTH', 'CURB', 'UCURB', 'PTOWB']
        for sp in ['H', 'D']:
            if 'BDENS_'+sp in self.file.variables.keys():   
                _keys = ['n_'+sp, 'orbloss_'+sp, 'cx1_'+sp, 'cx2_'+sp, \
                         'pi_'+sp, 'pe_'+sp, 'th_'+sp]
                _varnames = ['BDENS_'+sp,'BPLIM_'+sp, 'BPCXX_'+sp, 'BPCXI_'+sp,\
                             'PBI_'+sp, 'PBE_'+sp, 'PBTH_'+sp]
        for el in range(len(_keys)):
            keys.append(_keys[el])
            varnames.append(_varnames[el])
        self.nb_FKP_names, self.nb_FKP_vars = tu._fill_dict(self.file, keys, varnames)
        
        # Converting to metric system
        for k in keys:
            if k=='nbcd' or k=='unbcd':                
                self.nb_FKP_vars[k] *= 1e4
            
            elif k[0]=='n' or k[0:2]=='pi' or k[0:2]=='pe':
                self.nb_FKP_vars[k] *= 1e6
                
    def _calculate_all_nb(self):
        """first computations
        
        Hidden function to initialize and calculate variables
        i.e. vol-averaged temperatures and densities, ecc.

        Parameters:
            None
        Attributes:
            None
        Note:

        """
        self._calculate_scalars()
        self._calculate_eff()
        self._calculate_Ec()
        self._calculate_gi()
#        self._calculate_tau()
   
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
        A=2
        Aj = [2, self.imp_vars['A'].mean()]
        Zj = [1, self.imp_vars['Z'].mean()]
        nj = [self.kin_vars['ni'], self.imp_vars['nimp']]
        ne = self.kin_vars['ne']; Te = self.kin_vars['te']
        self.Ec = np.zeros(np.shape(self.rho) , dtype=float)
        self.ec_mean = np.zeros(self.nt, dtype=float)
        
        for it in range(self.nt):
            sum_facts = [nj[i][it,:]*(Zj[i]**2.)/Aj[i] for i in range(len(Aj))]
            summ = np.sum(np.asarray(sum_facts), axis=0)
            second_factor = (A**1.5/ne[it,:]*summ)**0.66
            Ec_prof = 14.8*Te[it,:]*second_factor
            self.Ec[it,:] = Ec_prof
            self.ec_mean[it] = np.trapz(Ec_prof, self.rho[it, :])
            
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
        self.taus = np.zeros((self.nt, len(self.rho[0,:])) , dtype=float)
        self.taus_mean = np.zeros(self.nt, dtype=float)
        
        if np.mean(ne)/1e18 >1:
            ne = ne*1e-6
        else:
            ne=ne

        for it in range(self.nt):
            ts=6.27e8*A[0]*np.divide(np.power(Te[it,:],1.5),(1.*lnlambda*ne[it,:]))
            taus=ts/3.*np.log((1+(E0_arr[it]/Ec[it,:])**1.5))           
            self.taus[it,:] = taus
            self.taus_mean[it] = np.dot(taus, self.dvol[it,:])/np.sum(self.dvol[it,:])

        
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
            self.gi_vavg[time_ind] = np.trapz(self.gi[time_ind,:], self.rho[time_ind, :])   
        self.gi_mean = np.mean(self.gi_vavg[self.inj_index])
        
    def _calculate_scalars(self):
        """ Scalar calculations

        Calculates scalar variables (integrals of profiles)

        Parameters:
            None
        Attributes:
            None    
        Note:
        """

        self.pe = np.zeros(self.nt, dtype=float)
        self.pi = np.zeros(self.nt, dtype=float)
        self.pth = np.zeros(self.nt, dtype=float)
        #difference between H and D beam
        self.Hpe = np.zeros(self.nt, dtype=float)
        self.Hpi = np.zeros(self.nt, dtype=float)
        self.Hpth = np.zeros(self.nt, dtype=float)
        self.Dpe = np.zeros(self.nt, dtype=float)
        self.Dpi = np.zeros(self.nt, dtype=float)
        self.Dpth = np.zeros(self.nt, dtype=float)

        self.nbcd = np.zeros(self.nt, dtype=float)
        self.unbcd = np.zeros(self.nt, dtype=float)
        self.jboot = np.zeros(self.nt, dtype=float)
        self.shield = np.zeros(self.nt, dtype=float)        

        for i in range(self.nt):
            self.pe[i] = np.dot(self.nb_FKP_vars['pe'][i,:], self.dvol[i,:])
            self.pi[i] = np.dot(self.nb_FKP_vars['pi'][i,:], self.dvol[i,:])
            self.pth[i] = np.dot(self.nb_FKP_vars['th'][i,:], self.dvol[i,:])
                
            self.nbcd[i] = np.dot(self.nb_FKP_vars['nbcd'][i,:], self.darea[i,:])
            self.unbcd[i] = np.dot(self.nb_FKP_vars['unbcd'][i,:], self.darea[i,:])           
            self.shield[i] = self.nbcd[i]/self.unbcd[i]           

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
        print("Power not counted anywhere:", pbal)
        print("Pbalance/Pin ", pbal/np.mean(self.nb_in_vars['P'])*1e6)
        print("psh+pol ", np.mean(self.pol+self.psh)*1e-6)
        print("pe+pi+pth ", np.mean(self.pe+self.pi+self.pth)*1e-6)
        print("pcx " , np.mean(self.pcxprof)*1e-6)


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
            self.nbcd
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


    def plot_deposition(self, time=[0]):
        """
        
        """
        plot_nb.plot_deposition(self, time)

    def plot_deposition_1d(self, time=[0], axp=0, axf=0, ls='-'):
        """
        """
        plot_nb.plot_deposition_1d(self, time, axp, axf, ls)

    def plot_deposition_prof(self, time=[0], axp=0, axf=0, ls='-'):
        """
        """
        plot_nb.plot_deposition_1d(self, time, axp, axf, ls)        
    
    def plot_radar_plost(self):
        """
        """
        self._calculate_scalars()
        N=3
        values = [self.psh_mean/self.pin_mean*100., self.pol_mean/self.pin_mean*100., \
                    self.pcx_mean/self.pin_mean*100.]

        categories = ['ST', 'OL', 'CX']
        au._radar_plot(N, values, categories,  title='Losses')

    def plot_NBCD(self, time=[0], axp=0, axf=0, ls='-'):
        """
        """
        plot_nb.plot_NBCD(self, time)        
        
        
    def plot_radar_deposition(self):
        """
        radar plot of absorbed power, NBCD, nf/ne, %powers to ions, efficiency, momentum        
        """
        N=4
        
        values = [self.pabs_mean/self.pin_mean*100., self.nbcd_mean*1e-3, \
                    self.nf_mean/self.ne_mean*100., self.gi_mean*100.]#, \
                    #self.eff_mean ]
        categories = [r'$P_{ABS}$', r'$I_{CD}$', r'$n_f/n_e$ [%]', r'$G_i$']#, r'$\eta$']
        au._radar_plot(N, values, categories, title='Deposition')
        
    def _init_pbalance_vars(self):
        """initialize vars for pbalance
        
        Parameters:
        Attributes:
            None    
        Note:        
        
        """
        x = self.t
        y = np.array([(self.pi+self.pe+self.pth)*1e-6, self.pcx*1e-6, self.pol*1e-6, self.psh*1e-6])
        yH = np.array([(self.Hpi+self.Hpe+self.Hpth)*1e-6, self.pcx*1e-6, self.pol*1e-6, self.psh*1e-6])
        pin = self.nb_in_vars['P']
        yfrac = y/pin*1e6
        return x,y,yH,yfrac
    
    def plot_powerbalance(self, time=[0], ls='-'):
        """ plot pow balance

        Plots the power balance with cumulative plot
       
        Parameters:
            time (arr) : array with the times   where to plot the lines
        Attributes:
            None    
        Note:
        """   
        plot_nb.plot_powerbalance(self, time, ls)
