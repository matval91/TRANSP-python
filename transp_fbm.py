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
from ascot_utils import common_style, limit_labels, _cumulative_plot, _plot_2d

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
        self.dict_dim['R'] *= 0.01
        self.dict_dim['z'] *= 0.01
        self.dict_dim_MCgrid = dict.copy(self.dict_dim)

        self.fdist_MCgrid = 0.5*self.infile.variables['F_D_NBI'][:]*1e6 #converting to m and normalizing over angle
        
        # This fdist_notnorm is the same as ascot distributions, this is useful because
        # I can copy the methods in this way       
        self.fdist_notnorm = np.zeros((self.shape_dim['E'], self.shape_dim['pitch'], \
                    self.nbins, self.nbins), dtype=float)        
        self._RZgrid()
        self.norm=1.
        self._computenorm()    
        self.fdist_norm = self.fdist_notnorm/self.norm
        self._readwall()

        self.rsurf, self.zsurf = self.infile.variables['RSURF'][:], self.infile.variables['ZSURF'][:]
        self.rsurf *= 0.01
        self.zsurf *= 0.01
        
        
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
        self.shot = self.runid[0:5]
        
        
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
        dist_toint = self.fdist_notnorm[:,:,:,:]

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #pitch,E, normalized           


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


    def _integrate_Ep(self):
        """
        Hidden method to integrate over (E,p)
        """
        dist_toint = self.fdist_notnorm/self.norm
            
        int_E = np.trapz(dist_toint, self.dict_dim['E'], axis=0)
        self.f_Ep_int = np.trapz(int_E, self.dict_dim['pitch'], axis=0)

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
        plot 2D (pitch, energy, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        x,y = self.dict_dim['pitch'], self.dict_dim['E']
        title = self.runid + ' ' + str(self.time)
        self._plot_2d(x,y, self.f_space_int, r'$\xi$', r'E [eV]', title)

    def plot_space(self):
        """
        """
        try:
            self.f_Ep_int.mean()
        except:
            self._integrate_Ep()
        x,y = self.dict_dim['R'], self.dict_dim['z']
        title = self.runid + ' ' + str(self.time)
        self._plot_2d(x,y, self.f_Ep_int.T, r'R [m]', r'z [m]', title, \
                    wall=[self.R_w, self.z_w], surf=[self.rsurf, self.zsurf])


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
        
    def _plot_2d(self,x,y,z, xlabel, ylabel, title, **kwargs):
        """
        """        
        fig = plt.figure()
        ax = fig.gca()
        if "wall" in kwargs:
            r,zt = kwargs['wall'][0], kwargs['wall'][1]
            ax.plot(r,zt, 'k', lw=2.5)
        if "surf" in kwargs:
            r,zt = kwargs['surf'][0], kwargs['surf'][1]
            ind=np.linspace(0, np.shape(r)[0]-1, 5)            
            for i in ind:
                plt.plot(r[i,:], zt[i,:], 'k', lw=1.1)
        CS = ax.contourf(x,y,z, 20,  cmap=my_cmap)
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
            print fname
            tmp = fbm(fname)
            self.timeslices[i] = tmp

        self._integratespace()
        self._integrateEp()
        self.R_w, self.z_w = self.timeslices[0].R_w, self.timeslices[0].z_w
        
        
    def _integratespace(self):
        """
        """
        for i in self.timeslices:
            i._integrate_space()
    
    def _integrateEp(self):
        """
        """
        for i in self.timeslices:
            i._integrate_Ep()
    
    def plot_Epframes(self):
        """
        """
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)        
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
                    self.timeslices[i].dict_dim['E']*1e-3, self.timeslices[i].f_space_int, \
                    20, cmap=my_cmap, vmin=0, vmax=vmax)
            #plt.colorbar(CB)        
            axarr[ind_row, ind_col].set_xlabel(r'$\xi$')
            axarr[ind_row, ind_col].set_ylabel(r'E [keV]')
            axarr[ind_row, ind_col].set_xlim([-1., 1.])
            axarr[ind_row, ind_col].set_ylim([0, 30])
            axarr[ind_row, ind_col].set_title(str(self.timeslices[i].time))
        
        f.tight_layout()
        f.subplots_adjust(right=0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.colorbar(CB, cax=cbar_ax)        
        f.text(0.03, 0.02, self.runid)
        plt.show()


    def plot_spaceframes(self):
        """
        """
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        n_cols = 4
        n_rows = self.nslices/n_cols
        if self.nslices%n_cols!=0:
            n_rows=n_rows+1
        vmax = np.max([i.f_Ep_int for i in self.timeslices])
        f, axarr = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
        xlim, ylim = [0.5, 1.2], [-0.8, 0.8]
        f.suptitle(r'Normalized fast ion distribution function', fontsize=16)
        for i in range(self.nslices):
            ind_row = i/n_cols
            ind_col = i%n_cols
            CB=axarr[ind_row, ind_col].contourf(self.timeslices[i].dict_dim['R'], \
                    self.timeslices[i].dict_dim['z'], self.timeslices[i].f_Ep_int.T, \
                    20, cmap=my_cmap, vmin=0, vmax=vmax)
            axarr[ind_row,ind_col].plot(self.R_w, self.z_w, 'k', lw=3.)
            #plt.colorbar(CB)        
            axarr[-1, ind_col].set_xlabel(r'R [m]')
            axarr[ind_row, 0].set_ylabel(r'z [m]')
            start, end = axarr[ind_row,0].get_ylim()
            axarr[ind_row,0].yaxis.set_ticks(np.arange(start, end, 0.4)) 
            start, end = axarr[-1, ind_col].get_xlim()
            axarr[-1, ind_col].xaxis.set_ticks(np.arange(start, end, 0.3)) 
            
            axarr[ind_row, ind_col].set_xlim(xlim)
            axarr[ind_row, ind_col].set_ylim(ylim)
            axarr[ind_row, ind_col].set_title(str(self.timeslices[i].time))
            r, zt = self.timeslices[i].rsurf, self.timeslices[i].zsurf            
            ind=np.linspace(0, np.shape(r)[0]-1, 5)            
            for j in ind:
                axarr[ind_row,ind_col].plot(r[j,:], zt[j,:], 'k', lw=1.1)
       
        f.tight_layout()
        f.subplots_adjust(right=0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.colorbar(CB, cax=cbar_ax)        
        f.subplots_adjust(top=0.88)
        f.text(0.03, 0.02, self.runid)
        plt.show()


#    def plot_frames(self):
#        """
#        """
#        
#        n_cols = 3
#        n_rows = self.nslices/3
#        if self.nslices%3!=0:
#            n_rows=n_rows+1
#        vmax = np.max([i.f_Ep_int for i in self.timeslices])
#        f, axarr = plt.subplots(n_cols, n_rows, sharex=True, sharey=True)        
#        for i in range(self.nslices):
#            ind_row = i/n_rows
#            ind_col = i%n_rows
#            CB=axarr[ind_row, ind_col].contourf(self.timeslices[i].dict_dim['R'], \
#                    self.timeslices[i].dict_dim['z'], self.timeslices[i].f_Ep_int, \
#                    20, cmap=my_cmap, vmin=0, vmax=vmax)
#            #plt.colorbar(CB)        
#            axarr[ind_row, ind_col].set_xlabel(r'R [m]')
#            axarr[ind_row, ind_col].set_ylabel(r'z [m]')
#            axarr[ind_row, ind_col].set_xlim([-0.5, 1.6])
#            axarr[ind_row, ind_col].set_ylim([-0.8, 0.8])
#            axarr[ind_row, ind_col].set_title(str(self.timeslices[i].time))
#            axarr[ind_row, ind_col].plot(self.R_w, self.z_w, 'k', lw=2.5)
#        f.tight_layout()
#        f.subplots_adjust(right=0.8)
#        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
#        f.colorbar(CB, cax=cbar_ax)        
#        
#        plt.show()

       
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
            title = "Shot "+str(el.shot)+'| t='+str(el.time)+'s'

            xlim, ylim = [-1., 1.], [0, 30000]
            _save_singleframe(x, y, z, xlab, ylab, title, xlim, ylim, vmax, str(i))
        os.chdir(olddir)


    def make_gif_space(self):
        """
        http://superfluoussextant.com/making-gifs-with-python.html
        """
        cwd = os.getcwd()
        tmpdir = cwd+'tmp_gifcreator/'
        self._make_space_singlefigures(tmpdir)
        os.chdir(tmpdir)
        
        gif_name = cwd+'space'+self.runid
        file_list = glob.glob('*.png') # Get all the pngs in the current directory
        #list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0])) # Sort the images by #, this may need to be tweaked for your use case

        with open('image_list.txt', 'w') as file:
            for item in file_list:
                file.write("%s\n" % item)
            
        os.system('convert -delay 100 @image_list.txt {}.gif'.format(gif_name)) # On windows convert is 'magick'
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        
    def _make_space_singlefigures(self, tmpdir):
        """
        """
        olddir = os.getcwd()
        os.makedirs(tmpdir); os.chdir(tmpdir)

        vmax = np.max([i.f_Ep_int for i in self.timeslices])
        for i, el in enumerate(self.timeslices):
            x,y = el.dict_dim['R'], el.dict_dim['z']
            z = el.f_Ep_int.T
            xlab, ylab = r'R [m]', r'Z [m]'
            title = "Shot "+str(el.shot)+'| t='+str(el.time)+'s'

            xlim, ylim = [0.5, 1.2], [-0.8, 0.8]
            _save_singleframe(x, y, z, xlab, ylab, title, xlim, ylim, vmax, str(i) ,\
                            wall_rz = [el.R_w, el.z_w], surf=[el.rsurf, el.zsurf])
        os.chdir(olddir)      
        
def _save_singleframe(x, y, z, xlab, ylab, title, xlim, ylim, vmax, savename, **kwargs):
    """
    """
    interactive(False) 
    levels = np.linspace(0,vmax,30)
    f = plt.figure()
    ax = f.add_subplot(111)        
    CB=ax.contourf(x,y,z, 20, cmap=my_cmap, vmin=0, vmax=vmax, levels=levels)
    plt.colorbar(CB)   
    if 'wall_xy' in kwargs:
        theta=np.linspace(0, math.pi*2., num=100)
        wall = kwargs['wall_xy']
        ax.plot(np.min(wall[0])*np.cos(theta), np.min(wall[0])*np.sin(theta), 'k', lw=2.5)
        ax.plot(np.max(wall[0])*np.cos(theta), np.max(wall[0])*np.sin(theta), 'k', lw=2.5)
    if 'wall_rz' in kwargs:
        wall = kwargs['wall_rz']
        ax.plot(wall[0], wall[1], 'k', lw=2.5)
        ax.axis('equal')

    if 'surf' in kwargs:
        r,zt = kwargs['surf'][0], kwargs['surf'][1]
        ind=np.linspace(0, np.shape(r)[0]-1, 5)            
        for i in ind:
            ax.plot(r[i,:], zt[i,:], 'k', lw=1.1)

    ax.grid('on')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title)
    plt.savefig(savename+'.png', bbox_inches='tight')
    interactive(True)
    
def _radar_plot(N, values, categories, title):
    """
    Freely revisited from 
    https://python-graph-gallery.com/390-basic-radar-chart/
    """    
    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values += values[:1]
    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * math.pi for n in range(N)]
    angles += angles[:1]
    angles = np.array(angles)
     
    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)
     
    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='black', size=30)
     
    # Draw ylabels
    ax.set_rlabel_position(45)
    #ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
#    plt.yticks([10,20,30], ["10","20","30"], color="grey", size=7)
#    plt.ylim(0,40)

    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid')
     
    # Fill area
    ax.fill(angles, values, 'b', alpha=0.1)
    
    ax.set_title(title, va='bottom')
