#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:43:58 2018

@author: vallar
"""
from __future__ import print_function
import netCDF4 as nc
import numpy as np
from matplotlib import interactive
import scipy.interpolate as interp
import glob, os, shutil
from utils.plot_utils import common_style, _cumulative_plot, _plot_2d
common_style()

class absorption:
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
    def __init__(self, inf_name, fname_surf=''):
        """
        Gathers variables
        """
        self.infile_n = inf_name
        self.file = nc.Dataset(inf_name)
        self.fname_surf = fname_surf
        self.shot = inf_name[0:5]
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
        self.data_i['phi']=self.data_i['phi']*np.pi/180.

        #convert the distances from cm to m
        for i in ['R','z', 'Rgc', 'zgc']:
            self.data_i[i] *= 0.01
            
        self.data_i['X'] = np.multiply(self.data_i['R'],np.cos(self.data_i['phi']))
        self.data_i['Y'] = np.multiply(self.data_i['R'],np.sin(self.data_i['phi']))

        self.time = self.file.variables['bs_time_D_MCBEAM'][:].mean().round(3)
        self.read_surf()
        self._readwall()
            
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


    def read_surf(self):
        """
        """
        if self.fname_surf!='':
            self.infile_surf = nc.Dataset(self.fname_surf)
            rsurf, zsurf = self.infile_surf.variables['RSURF'][:], self.infile_surf.variables['ZSURF'][:]
            _rho =  self.infile_surf.variables['XSURF'][:]
            rsurf *= 0.01; zsurf *= 0.01
            self.surf= np.array([rsurf, zsurf])
            rho = np.repeat(_rho, (np.shape(rsurf)[1]))   
            grid_x, grid_y = np.mgrid[np.min(rsurf):np.max(rsurf):100j, np.min(zsurf):np.max(zsurf):200j]
            self.RZsurf = interp.griddata(np.array([rsurf.flatten(), zsurf.flatten()]).T, rho, (grid_x, grid_y), fill_value=1) 
            self.RZsurf_param = interp.interp2d(grid_x, grid_y, self.RZsurf)
        else:
            self.surf=0;
        
    def plot_RZ(self):
        """
        Method to plot R vs z of the ionised particles, useful mostly with bbnbi
        """
        x=self.data_i['R']
        y=self.data_i['z']
        wallrz= [self.R_w, self.z_w]

        xlab = 'R [m]'
        ylab = 'z [m]' 
        _plot_2d(x, y, xlabel=xlab, ylabel=ylab,\
                dist=0, title='t={:.2f} s'.format(self.time), wallxy=0, wallrz=wallrz, surf=self.surf, R0=0, ax=0, \
                scatter=0, hist=1, xlim=[0.6, 1.1], ylim=0, fname='', cblabel='', lastpoint=1)        
        
    def plot_XYpart(self):
        """
        Method to plot XY of ionisation, without difference between the beams
        """
        x=self.data_i['X']
        y=self.data_i['Y']
        wallrz= [self.R_w, self.z_w]
        
        _plot_2d(x,y, xlabel=r'X (m)', ylabel=r'Y (m)',\
           dist=0, title='t={:.2f} s'.format(self.time), wallxy=wallrz, wallrz=0, surf=0, R0=0.88, ax=0, \
           scatter=1, hist=0, xlim=0, ylim=0, fname='', cblabel='', lastpoint=1)       
    
    def plot_Epitch(self):
        """
        Method to plot E pitch of ionisation
        """
        y=self.data_i['E']
        x=self.data_i['pitch']

        _plot_2d(x,y, xlabel=r'$\xi$', ylabel=r'E (keV)',\
                dist=0, title='t={:.2f} s'.format(self.time), wallxy=0, wallrz=0, surf=0, R0=0, ax=0, \
                scatter=0, hist=1, xlim=0, ylim=0, fname='', cblabel='', lastpoint=1)            
        
    def _readwall(self):
        """
        Hidden method to read the wall
        """
        in_w_fname='/home/vallar/TCV_wall/TCV_vessel_coord.dat'
        try:
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=0)
        except OSError:
            print('No wall available')
            wall=np.zeros((2,1))
        self.R_w = np.array(wall[0,:])
        self.z_w = np.array(wall[1,:])
        
class absorption_time:
    """
    Class as fbm including time dependence
    """
    def __init__(self, runid):
        self.runid = runid
        file_list = glob.glob('*_birth.cdf*') # Get all the cdf in the current directory
        self.nslices=len(file_list)        
        self.timeslices = np.array([absorption for _ in range(self.nslices)])
        print("Number of timeslices: ", str(self.nslices))
        
        
        for i in range(self.nslices):
            fname = runid+'_birth.cdf'+str(i+1)
            self.nslices += 1
            print(fname)
            tmp = absorption(fname)
            self.timeslices[i] = tmp
        self.R_w, self.z_w = \
        self.timeslices[0].R_w,self.timeslices[0].z_w

    def make_gif_XY(self):
        """
        """
        self._make_gif_coords('X', 'Y', 'X [m]', 'Y [m]', wall=[self.R_w, self.z_w])

    def _make_gif_coords(self, xname, yname, xlab, ylab, wall):
        """
        http://superfluoussextant.com/making-gifs-with-python.html
        """
        cwd = os.getcwd()
        tmpdir = cwd+'tmp_gifcreator/'      
        self._make_coords_singlefigures(tmpdir, xname, yname, xlab, ylab, wall)
        os.chdir(tmpdir)
        
        gif_name = cwd+xname+yname+'_ioniz_'+self.runid
        file_list = glob.glob('*.png') # Get all the pngs in the current directory
        #list.sort(file_list, key=lambda x: int(x.split('_')[1].split('.png')[0])) # Sort the images by #, this may need to be tweaked for your use case

        with open('image_list.txt', 'w') as file:
            for item in file_list:
                file.write("%s\n" % item)
            
        os.system('convert -delay 100 @image_list.txt {}.gif'.format(gif_name)) # On windows convert is 'magick'
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        
    def _make_coords_singlefigures(self, tmpdir, xname, yname, xlab, ylab, wall):
        """
        """
        olddir = os.getcwd()
        os.makedirs(tmpdir); os.chdir(tmpdir)
    
        vmax = 15
        for i, el in enumerate(self.timeslices):
            x = np.linspace(np.min(el.data_i[xname]), np.max(el.data_i[xname]), 100)
            y = np.linspace(np.min(el.data_i[yname]), np.max(el.data_i[yname]), 100)
            title = "Shot "+str(el.shot)+'| t='+str(el.time)+'s'
            interactive(False) 
            z, xedges, yedges = np.histogram2d(el.data_i[xname], el.data_i[yname],\
                            bins=100, normed=True)     
            z=z.T
            xlim, ylim = [-max(self.R_w)*1.1, max(self.R_w)*1.1], [-max(self.R_w)*1.1, max(self.R_w)*1.1]
            interactive(True) 
            _save_singleframe(x, y, z, xlab, ylab, title, xlim, ylim, vmax, str(i), \
                wall_xy = [self.R_w])
        os.chdir(olddir)


class particles:
    """
    Handles the particles position
     u'nstat_track_xja_D_NBI',
     
     u'track_phiin_D_NBI',
     u'track_phiout_D_NBI',

     u'track_rin_D_NBI',
     u'track_rout_D_NBI',

     u'track_zin_D_NBI',
     u'track_zout_D_NBI',

     u'track_einj_D_NBI',
     u'track_sinj_D_NBI', #weight

     u'track_xl_D_NBI',
     u'track_pdep_D_NBI',
     u'track_vpv_dep_D_NBI',
     u'track_dr_flr_D_NBI',
     u'track_dz_flr_D_NBI',
     u'track_de_flr_D_NBI',
     u'track_dvpv_flr_D_NBI',

    
    """
    def __init__(self, infile_n):
        self.infile_n = infile_n
        self.infile = nc.Dataset(infile_n)
        self._read_info()
        self.dictin = {'R':'track_rin_D_NBI', 'phi':'track_phiin_D_NBI', \
                    'z':'track_zin_D_NBI', 'pitch':'A', 't':'TIMELST'}
        
        
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
        self._readwall()
