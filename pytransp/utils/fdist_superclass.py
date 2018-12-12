#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
import numpy as np
import math
import utils.plot_fdist as plot_fdist

class fbm_superclass:
    """
 
    """
    def __init__(self, infile_n):
        """       
        """
        self.fname = infile_n
        self._read()
        
    def _computenorm(self):
        """
        calculates the norm of the function
        """
        if "pitch" in self.dict_dim.keys():
            try:
                self.f_spacep_int()
            except:
                self._integrate_spacep()
                
            self.norm = np.trapz(self.f_spacep_int, self.dict_dim['E'])
            self.norm=1.
            #print "NORM = ", self.norm

    def _readwall(self):
        """
        Hidden method to read the wall
        """
        in_w_fname='/home/vallar/TCV/TCV_vessel_coord.dat'
        wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=0)

        self.R_w = np.array(wall[0,:])
        self.z_w = np.array(wall[1,:])        
       
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
            
    def _integrate_space_rz(self):
        """
        Function to integrate over (R,z)
        """
        dist_toint = self.fdist_notnorm[:,:,:,:]

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #pitch,E, normalized           
             
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


    def _integrate_Ep(self):
        """
        Hidden method to integrate over (E,p)
        """
        dist_toint = self.fdist_notnorm/self.norm
            
        int_E = np.trapz(dist_toint, self.dict_dim['E'], axis=0)
        self.f_Ep_int = np.trapz(int_E, self.dict_dim['pitch'], axis=0)
        
    def _integrate_spacemu(self):
        """
        hidden method to integrate over (space,mu)
        """
        self.f_spacemu_int = self._integrate_spacex('mu')        
        
        
        
    def plot_spaceE(self, ax=0, color='', ls='-'):
        """
        plot 1D (pitch, int_space (int_E (fdist)))
        """
        try:
            self.f_spaceE_int.mean()
        except:
            self._integrate_spaceE()    
        plot_fdist.plot_spaceE(self, ax, color, ls)
        
    def plot_spacep(self, ax=0, color='', ls='-'):
        """
        plot 1D (energy, int_space (int_pitch (fdist)))
        """
        try:
            self.f_spacep_int.mean()
        except:
            self._integrate_spacep()
            
        plot_fdist.plot_spaceE(self, ax, color, ls)

    def plot_Epitch(self, ax=0):
        """
        plot 2D (pitch, energy, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        plot_fdist.plot_Epitch(self, ax)

    def plot_space(self, ax=0):
        """
        """
        try:
            self.f_Ep_int.mean()
        except:
            self._integrate_Ep()
        plot_fdist.plot_space(self, ax)
