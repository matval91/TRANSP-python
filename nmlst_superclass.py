#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:40:49 2018

@author: vallar
"""
import scipy.interpolate as interp
import numpy as np
import tcv
import os
import read_sig


class NMLST:
    """
    Class defining all the variables needed to be set in the namelist
    This is the list of values needed:
        $SHOT, $TBEG, $TEND, $DT
    
    """
    def __init__(self, indict):
        self.indict = indict
        self.shot = indict['shot']
        self.run  = indict['run']
        self.path = indict['path']

        self.tbeg = indict['tbeg']
        self.tend = indict['tend']
        self.dt   = indict['dt']
        self.t = np.linspace(self.tbeg, self.tend, self.dt, dtype=float)
#    replace_dict = {'$NPTCLS':self.nptcls,\
#                    '$SHOT': self.shot,'$timestep':self.dt,\
#                    '$tinit': self.tbeg, '$ftime':self.tend,\
#                    '$NLJCCW': int(self.jsgn), '$NLBCCW': int(self.bsgn),\
#                    '$EINJ': self.einj, '$ZBEAM':self.zbeam, '$ZEFF':self.zeff}    
    def _read_nmlst(self, fname):
        """
        
        """
        with open(fname) as f_nml:
            lines_nml_tmp=f_nml.readlines()
        
        
        lines_nml=''.join(lines_nml_tmp)
        lines_nml=lines_nml.replace('\r','')        
        return lines_nml

    def _interp1d(self,x,y):
        """
        interpolates 1D data to the time requested
        """
        dummy = interp.InterpolatedUnivariateSpline(x, \
                                                    y, ext=0)
        newy = dummy(self.t)
        return newy
    
    def merge(self):
        """
        Merges with existing namelist
        """
        if not os.path.isfile(self.std_fname):
            print "No std namelist existent ", self.std_fname
            return
        
        lines_std = self._read_nmlst(self.std_fname)
        lines_add = self._read_nmlst(self.out_fname)
        lines_out = lines_std+lines_add
        self._backup_namelist(self.std_fname)
        self._backup_namelist(self.out_fname)
        outfile=open(self.std_fname,'w')
        outfile.write(lines_out)
        outfile.close()
        print "Merged ", self.std_fname, ' and ', self.out_fname, ' in ', self.std_fname
        
    def write(self):
        """
        """
        print "Writing namelist"
        lines_nml = self._read_nmlst(self.template_fname)
            
        self.lines_out=lines_nml
        for i,key in enumerate(self.rep_dict):
            print key, self.rep_dict[key]
            if np.isscalar(self.rep_dict[key]):
                self.lines_out=self.lines_out.replace\
                         (\
                          key,str(self.rep_dict[key])\
                          )
            else:
                #This applies for arrays
                self.lines_out=self.lines_out.replace\
                         (\
                          key, np.array2string(self.rep_dict[key], separator=',')[1:-1]\
                          )
        self._backup_namelist(self.out_fname)
        outfile=open(self.out_fname,'w')
        outfile.write(self.lines_out)
        outfile.close()
    
    def _backup_namelist(self, fname):
        """
        Does backup of a namelist adding ~
        """
        if os.path.isfile(fname):
            os.rename(self.out_fname, self.out_fname+'~')
            print "\n Copy ", self.out_fname, " to ", self.out_fname+'~'            
    
    def _JBsign(self):
        """
        Reads sign of I and B using read_sig class
        """
        sig = read_sig.DWR1(self.indict)
        self.jsgn_b = 'T'
        self.bsgn_b = 'T'
        self.jsgn, self.bsgn = sig._readJB()
        self.jb = self.jsgn*self.bsgn
        if self.jsgn<0:
            self.jsgn_b = 'F'
        if self.bsgn<0:
            self.bsgn_b = 'F'        
            
class basic_nmlst(NMLST):
    """
    Standard namelist, with the following values needed:
        $SHOT, $TBEG, $TEND, $DT, $BSGN, $JSGN       
    """
    def __init__(self, indict):
        NMLST.__init__(self, indict)
        self.rep_dict = {'$SHOT':indict['shot'], '$TBEG':indict['tbeg'],\
                         '$TEND':indict['tend'], '$DT':indict['dt'],\
                         '$JSGN':'', '$BSGN':''}
        self.template_fname = '/home/vallar/matlab_TRANSP/namelist_template_py.DAT'
        self.out_fname = self.path+str(self.shot)+str(self.run)+'TR.DAT'
        self._JBsign()
        self.rep_dict['$JSGN'] = self.jsgn_b
        self.rep_dict['$BSGN'] = self.bsgn_b
        
        
        
class NBI_nmlst(NMLST):
    """
    NBI namelist, with the following values needed:
    the letters D and H refer to the species injected (NBH:D, diagnostic beam:H)
        $FBM_outtim 
        $NBPTCLS
        $ZBEAM_D:
        $EINJ_DBEAM
        $FFULL_DBEAM
        $FALF_DBEAM, 
        $EINJ_HBEAM
        $FFULL_HBEAM
        $FHALF_HBEAM
    """    
    def __init__(self, indict, **kwargs):
        NMLST.__init__(self, indict)
        self.rep_dict = {'$FBM_outtim': np.array([0]),'$NBPTCLS':80000,'$ZELEV_D':0.,\
        '$EINJ_DBEAM':0., '$FFULL_DBEAM':0., '$FHALF_DBEAM':0., \
        '$EINJ_HBEAM':0., '$FFULL_HBEAM':0., '$FHALF_HBEAM':0.}  
        
        for key in kwargs:
            if key=='FBM_outtim':
                self.rep_dict['$'+key] = np.array(kwargs[key])                         
            else:
                self.rep_dict['$'+key] = kwargs[key]         

        self.template_fname = '/home/vallar/matlab_TRANSP/namelist_NBI_template.DAT'
        self.std_fname = self.path+str(self.shot)+str(self.run)+'TR.DAT'
        self.out_fname = self.std_fname+'_NBI'
        
        self.conn = tcv.shot(indict['shot'])
        self._fill_dict()
    def _fill_dict(self):
        """
        Reads data from tree
        """
        self.get_NB_heating_data()
        self.get_NB_diag_data()
        for key in self.rep_dict.keys():
            print key, self.rep_dict[key]
    
    def get_NB_heating_data(self):
        """
        Reads NB heating data
        """
        self._NB_heating_energy()
        self._NB_heating_fractions()
        
    def _NB_heating_energy(self):
        """
        reads energy of the beam
        """
        nbi = self.conn.tdi(r'\ATLAS::NBH.DATA.MAIN_ADC:DATA')
        tim = nbi.dim_0
        acc_volt = np.asarray(nbi[:,32]) # acc voltage in kV
        _ind = np.where(np.logical_and(\
                                      np.logical_and(tim>self.tbeg, tim<self.tend),
                                      acc_volt>0.2*max(acc_volt)
                                      )
                         )
        self.einj = np.mean(acc_volt[_ind])*1e3
        self.rep_dict['$EINJ_DBEAM'] = round(self.einj,2)
        
        
    def _NB_heating_fractions(self):
        """
        Heating beam energy fraction
        """
        try:
            self.einj.mean()
        except:
            self._NB_heating_energy()
        
        try:
            self.conn+1. #just to give an error
            # TO ADD SIGNALS TO READ
        except:
            if self.einj>20e3:
                self.rep_dict['$FFULL_DBEAM'] = 0.76
                self.rep_dict['$FHALF_DBEAM'] = 0.20
            else:
                self.rep_dict['$FFULL_DBEAM'] = 0.50
                self.rep_dict['$FHALF_DBEAM'] = 0.25
                
    def get_NB_diag_data(self):
        """
        Getting data for diagnostic beam (now set as 58823)
        """
        self.rep_dict['$EINJ_HBEAM']  = 49828.00
        self.rep_dict['$FFULL_HBEAM'] = 0.6378
        self.rep_dict['$FHALF_HBEAM'] = 0.0712
         
         
class EC_nmlst(NMLST):
    """
    EC namelist, with the following values needed:
        $NANTECH  ! number of launchers
        $XECECH ! MAJOR RADIUS OF SOURCE (CM)
        $ZECECH ! HEIGHT ABOVE MIDPLANE OF
        $FREQECH ! WAVE FREQUENCY (Hz)
        $PHAIECH ! AZIMUTHAL LAUNCH ANGLE (toroidal)
        $THETECH ! POLAR LAUNCH ANGLE (theta_pol)

        $RFMODECH 0.0, 0.0, 0.0, 0.0, 0.0, ! power fraction in o-mode (0 or 1)
        $EFFECH 1.0, 1.0, 1.0, 1.0, 1.0, ! effective current drive (fudge factor)
        $NRAYECH 30,30,30,30,30, ! number or rays
        $BSRATECH 30.0, 30.0, 30.0, 30.0, 30.0, ! Aspect ratio of beam (1 for circle)
        $BHALFECH 6.0, 6.0, 6.0, 6.0, 6.0, ! divergence angle (degree)
        $NDAMPECH 2, 2, 2, 2, 2, ! damping model (good for x2, take 3 for x3)

    """    
    def __init__(self, indict, **kwargs):
        NMLST.__init__(self, indict)
        self.rep_dict = {'$NANTECH':1,'$XECECH':[0.],'$ZECECH':[0.0],\
        '$FREQECH':[0.], '$PHAIECH':[0.], '$THETECH':[0.], \
        '$RFMODECH':[0.], '$EFFECH':[0.], '$NRAYECH':[0.],\
        '$BSRATECH':[0.], '$BHALFECH':[0.], '$NDAMPECH':[0.]}  
       
        self._TCV_ECparam()
        self.template_fname = '/home/vallar/matlab_TRANSP/namelist_EC_template.DAT'
        self.std_fname = self.path+str(self.shot)+str(self.run)+'TR.DAT'
        self.out_fname = self.std_fname+'_EC'
        
        sig_ec = read_sig.DWR(indict)
        sig_ec._read_ecrh()
        self.ind_gyro = sig_ec.ind_Gyro
        self.jsgn, self.bsgn = sig_ec._readJB()
        self.jb = self.jsgn*self.bsgn
        
        self.rep_dict['$NANTECH']  = len[self.ind_gyro]
        self.rep_dict['$FREQECH']  = self.frequency[self.ind_gyro]
        self.rep_dict['$RFMODECH'] = self.omode[self.ind_gyro]
        self.rep_dict['$EFFECH']   = np.full(len(self.ind_gyro), 1.0, dtype=float)
        self.rep_dict['$NRAYECH']  = self.raysbylauncher[self.ind_gyro] 
        self.rep_dict['$BHALFECH'] = np.full(len(self.ind_gyro), 6.0, dtype=float)
        self.rep_dict['$NDAMPECH'] = self.ndampech[self.ind_gyro]
        
        self.conn = tcv.shot(indict['shot'])
        self._read_EC_geom()

    def _read_EC_geom(self):
        """
        Read from tree the data for X2 and X3
        """         
#        time = self.tbeg+(self.tend-self.tbeg)/2.
#        
#        [r_phi0,z_phi0,r1,z1,rc,zc,thetphi_L,thetphi_tor,raysbylauncher,powers,fname,z_axis] = \ 
#        toray_raygeom_TCV(self.shot,time,self.frequency,\
#                          self.omode,self.raysbylauncher,self.liuqe_vers,\
#                          self.path,self.isave,self.ipowers_tree,self.jb)
        

        
    def _TCV_ECparam(self):
        """
        Data for EC antennas in TCV (01/2018)
        """
        self.frequency = [82.7e9, 82.7e9, 82.7e9, 82.7e9, 82.7e9, 82.7e9, \
                          118.e9, 118.e9, 118.e9]
        self.omode = [0., 0., 0., 0., 0., 0., 0., 0., 0.] # o mode fraction
        self.nrays=30;
        self.raysbylauncher=np.full(9, self.nrays, dtype=int);    
        self.ndampech = [2,2,2,2,2,2,3,3,3]
        self.liuqe_vers=1; #version 1! is that ok?
        self.isave=1; #no save file
        self.ipowers_tree= 1; #=> read powers from tree