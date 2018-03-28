
from __future__ import print_function
import numpy as np
#GUI
import Tkinter as tk
# for reading and writing ufiles
import ufiles
# for the definition of appropriate plotting style
import tksty
# for the appropriate computation of momentum of equilibrium
import mom2rz
# for properly handling the equilibrium
import eqtools
#from tooltip import createToolTip
import sys
import contourc
# this is used for proper description of
# used as a fast way to fit the data
import descur 
# interpolation for increasing resolution
# in (R, Z)
from scipy import interpolate, integrate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg

import ReadEQDSK
#matplotlib.use("TkAgg")

sty = tksty.TKSTY()
descu  = descur.DESCUR().descur_fit
if sys.platform =='darwin':
    sys.path.append('/Users/vianello/Documents/Fisica/Computing/pythonlib/signalprocessing')
else:
    sys.path.append('/home/vianello/pythonlib/signalprocessing')

class GEO:

    """
    ====
    GEO SUPER Class
    ====

    This GEO class defines the common function among the two inherited classes indict and eqdsk.
    Input:
    ---------------

    A dictionary containing all the needed information. The dictionary may
    be created through a GUI or directly given
    Dictionary is of the form
    {'shot': shot, 'tbeg': tbeg, 'tend':tend,'dt':time resolution, 
    'Nmom':number of fourier moments,'Mom type':4, 
    'n_rho': number of point in toroidal rho,'n_the': number of point in theta, 
    'run': run of the shot}
    
    e.g.
    inDict={'shot':53778, 'dt':0.02, 'tbeg':1.4, 'tend':1.4, 'n_rho':41, \
    'n_the':101, 'n_rho_eq':40, 'path': YOUR HOME, 'run':'A01'}

    The path where the ufiles are stored is automatically defined
    $HOME/tr_client/TCV/#SHOT/#RUN
    
    Methods:
    ---------------
        _storedict: stores the values passed as input in the dictionary
        sepx: return the separatrix (as recomputed accordingly to the moment
            representation in the (rho,theta) grid
        mag_axis: calculates the position of magnetic axis in time window
            chosen
        RZ_flux: Radial and Vertical grid of the toroidal flux
            equipotential surfaces
        calc_qprof: Return the Q profile on the toroidal grid
        calc_RBtor: Returns the poloidal flux function on toroidal grid
        store_u: stores rz, qprof and RBt in one turn
        store_rz: Store in a ufile R and Z (for the toroidal flux) as a
           function of (time, rhotor, theta)
        store_qprof: Store the safety factor profile in the rho toroidal grid
           in the appropriate ufile
        store_RBt: Store R.Bt in the rho toroidal grid
           in the appropriate ufile
        plot_input: plots 1D and 2D profiles for benchmark on data given

    Attributes:
    ---------------
        rz_grid: r, z grid used in the computation of the field
        rho: rho gridspec
        theta: theta gridspec
    """
    timeout = 1e-10
    tlbl = 'Time'.ljust(20) + 'Seconds'

    def __init__(self, indict):
        self._storedict(indict)
        print("\n")
        print("===================")
        print("Initialisation of Magnetics Done")
        print("===================")
        print("\n")
        
    def _storedict(self, indict):
        self.shot = indict['shot']
        self.n_rho = indict['n_rho_eq']
        self.n_the = indict['n_the']
        self.tbeg = indict['tbeg']
        self.tend = indict['tend']
        self.dt   = indict['dt']
        self.path = indict['path']
        self.mom_type = 4
        self.mom_order = 8
        
    def sepx(self):

        """ 
        It just gives back the representation of the separatrix
        from the previously calculated moment representation of the
        toroidal flux. Actually the (r,z)
        """

        # check if the r_plot and z_plot has not been inizialied otherwise
        # it calls the appropriate method
        try:
            self.moments.mean()
        except:
            self._eqflux()


        # TRICK TO PUT THE PLASMA TO 0
        #self.z_plot = self.z_plot - self.z_axis
        r_sp = self.r_plot[:, -1, :]
        z_sp = self.z_plot[:, -1, :]
        return r_sp, z_sp


    def mag_axis(self):
        """
        Returns magnetic axis position (R,Z) in time
        """
        # check of magnetic axis calculation
        try:
            self.r_axis.mean()
        except:
            self._mag_axis()

        return self.r_axis, self.z_axis


    def RZ_flux(self):
        """
        Return the R and Z coordinates of the toroidal flux iso-contour
        as calculated from moment representation, i.e as a function of
        (time, rho, theta)

        """
        # always check that eqflux has been computed
        try:
            self.moments.mean()
        except:
            self._eqflux()



        self.plot_input(self.r_plot, self.z_plot, 'RZ')
        print("Done RZ equilibrium")
        #return self.r_plot, self.z_plot

    def calc_qprof(self):
        """
        Return the Q profile on the toroidal grid as
        a function of time and rhotor
        Return the qprofile, toroidal flux maps
        """
        try:
            self.qprof.mean()
        except:
            self._calc_qprof()
                
        
        self.plot_input(self.qprof,'Q')
        print("Evaluated Q profile")
        #return self.qprof

    def calc_RBtor(self):
        """
        returns RBtor
        Method to compute R.Btor in as profile of (t,rho_tor)
        """

        try:
            self.RBT.mean()
        except:
            self._calc_RBtor()

        self.plot_input(self.RBT,'F')
        print("Evaluated F profile")

    def coord_change_torflux(self, torflux_in):
        """
        Applying formula
        phi(cocos_OUT) = sigma_b0*phi(cocos_IN)                        [TOROIDAL FLUX] 
        """
        try:
            self.sigma_ip+1
        except:
            self._coord_change()
        
        torflux_out = np.multiply(self.sigma_b0,torflux_in)
        return torflux_out
    
    def coord_change_F(self):
        """
        Applying formula 
        F  (cocos_OUT) = sigma_b0*F(cocos_IN)
        """     
        try:
            self.sigma_ip+1
        except:
            self._coord_change()
        tmp = np.multiply(self.RBT,self.sigma_b0)
        self.RBT = tmp

    def coord_change_q(self):
        """
        Applying formula 
        q  (cocos_OUT) = sigma_ip*sigma_b0*sigma_rtp*q(cocos_IN)
        """
        try:
            self.sigma_ip+1
        except:
            self._coord_change()
        fact = self.sigma_ip*self.sigma_b0*self.sigma_rtp
        tmp=np.multiply(fact, self.qprof)
        self.qprof = tmp


        
    def _coord_change(self):
        """
        Function to convert the COCOS
        LIUQE is in COCOS 17 while EFIT in COCOS 3
        i=COCOS 17, o=COCOS 3

        These formula come from appendix C of COCOS article
        psi(cocos_OUT) = sigma_ip*sigma_bp*(2pi)^(e_bp)*psi(cocos_IN)  [POLOIDAL FLUX]
        
        phi(cocos_OUT) = sigma_b0*phi(cocos_IN)                        [TOROIDAL FLUX]
        F  (cocos_OUT) = sigma_b0*F(cocos_IN)
        q  (cocos_OUT) = sigma_ip*sigma_b0*sigma_rtp*q(cocos_IN)
        """
        ############################################
        # FROM COCOS article, table 1
        e_bp_i = 1
        e_bp_o = 0
        self.e_bp = e_bp_o-e_bp_i

        sigma_bp_i = -1
        sigma_bp_o = -1
        self.sigma_bp = sigma_bp_i*sigma_bp_o

        #this is sigma R,phi,z: cylind
        sigma_rpz_i = +1
        sigma_rpz_o = +1
        self.sigma_rpz = sigma_rpz_i*sigma_rpz_o
        
        #this is sigma rho,theta,phi: poloid
        sigma_rtp_i = +1
        sigma_rtp_o = -1
        self.sigma_rtp = sigma_rtp_i*sigma_rtp_o
        ############################################

        ############################################
        # FROM APPENDIX C
        # a specific sign of Ip and Bp is not requested
        #self.sigma_ip = self.sigma_rpz
        #self.sigma_b0 = self.sigma_rpz
        
        ############################################

        ############################################
        # FROM APPENDIX C
        # a specific sign of Ip and Bp is requested (q and F must be >0)
        sign_ip_i = -1 # for shot 53778
        sign_ip_o = 1
        self.sigma_ip = sign_ip_i*sign_ip_o

        sign_b0_i = -1 #for shot 53778
        sign_b0_o = 1
        self.sigma_b0 = sign_b0_i*sign_b0_o
        
        ############################################
    
#    def store_sep_moments(self):
#        """
#        Store as u-file the appropriate moments definition of the
#        toroidal flux for the separatrix
#        """
#        try:
#            self.moments.mean()
#        except:
#            self._eqflux()
#
#        karr = 1 + np.arange(self.mom_type)
#        marr = 1 + np.arange(self.mom_order)
#        darr = self.moments[:, -1,  :, :]
#
#        klbl = 'MOMENT TYPE'
#        mlbl = 'MOMENT INDEX'
#        dlbl = 'FOURIER MOMENTS'.ljust(20) + 'M'
#
#        uf_d = {'pre': 'M', 'ext': 'MRY', 'shot': self.shot,
#                'grid': {'X': {'lbl': self.tlbl, 'arr': self.time_u},
#                         'Y': {'lbl': mlbl, 'arr': marr},
#                         'Z': {'lbl': klbl, 'arr': karr}},
#                'data': {'lbl': dlbl, 'arr': darr}}
#        ufiles.WU(uf_d, udir=self.path)


    def store_u(self):
        """
        Stores all the ufiles 
        """
        self.store_rz()
        self.store_qprof()
        self.store_RBt()

        
    def store_rz(self):
        """
        Store as u-file the reconstructed R, Z position of
        equi toroidal flux surfaces

        Parameters
        ----------
        None
        -------
        None

        """
        thlbl = 'THETA               RAD'
        rlbl  = 'MAJOR RADIUS        M'
        zlbl  = 'VERTICAL POSITION   M'
        tlbl = 'Time'.ljust(20) + 'Seconds   '
        rholbl = 'RHO_TOR'+' '*23
        
        uf_d = {'pre': 'R', 'ext': 'SURF', 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl, 'arr': self.time_u},
                         'Y': {'lbl': rholbl, 'arr': self.rho},
                         'Z': {'lbl': thlbl, 'arr': self.theta}},
                'data': {'lbl': rlbl, 'arr': self.r_plot}}
        ufiles.WU(uf_d, udir=self.path)

        uf_d = {'pre': 'Z', 'ext': 'SURF', 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl, 'arr': self.time_u},
                         'Y': {'lbl': rholbl, 'arr': self.rho},
                         'Z': {'lbl': thlbl, 'arr': self.theta}},
                'data': {'lbl': zlbl, 'arr': self.z_plot}}
        ufiles.WU(uf_d, udir=self.path)


    def store_mom(self):
        """
        Store the u-file of the moments of the isoflux surfaces

        Parameters
        ----------
        None
        -------
        None
        """        
        try:
            self.moments.mean()
        except:
            self._eqflux()
        # Now need to conver the 4D array of momentum [time, rho, order, type] to 3D array of momentum:
        # [time, rho, order*type].
        # order is the number of n*theta used (the, 2the, 3the,...,order*the), type is fixed to 4 since we associate
        # 2 periodic functions to the each of the 2 variables R,Z
        tmpmoments = np.transpose(self.moments, (0,1,3,2))
        self.moments_uf = np.reshape(tmpmoments, (self.nt, self.n_rho,\
                                                  self.mom_order*self.mom_type))
        self.moments_uf *= 100 #in centimeters
        tlbl   = 'Time'.ljust(20) + 'Seconds   '
        rholbl = ' RHOTOR'.ljust(30)
        ilbl   = ' MOMENT INDEX'.ljust(30)
        zlbl   = ' FOURIER MOMENTS'.ljust(20)+' CM'.ljust(10)
        
        uf_d = {'pre': 'C', 'ext': 'MMX', 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl, 'arr': self.time_u},
                         'Y': {'lbl': rholbl, 'arr': self.rho},
                         'Z': {'lbl': ilbl, 'arr': np.arange(self.mom_order*(self.mom_type))}},
                'data': {'lbl': zlbl, 'arr': self.moments_uf}}
        ufiles.WU(uf_d, udir=self.path)            
            
    def store_qprof(self):
        """
        Store the u-file of the reconstructured Q profile
        on the toroidal grid

        Parameters
        ----------
        None
        -------
        None
        """
        try:
            self.qprof.mean()
        except:
            self.calc_qprof()

        # now generate the appropriate dictionary for
        # saving the files

        tlbl = 'Time'.ljust(20) + 'Seconds   '
        rholbl = 'rho_tor'+' '*23
        lbl = 'Q profile'+' '*21
        uf_d = {'pre': 'Q', 'ext': 'SAF', 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl, 'arr': self.time_u},
                         'Y': {'lbl': rholbl, 'arr': self.rho}},
                'data': {'lbl': lbl, 'arr': self.qprof}}
        ufiles.WU(uf_d, udir=self.path)


    def store_RBt(self):
        """
        Store the u-file of the R.Bt profile
        on the toroidal grid

        Parameters
        ----------
        None
        -------
        None
        """
        try:
            self.RBT.mean()
        except:
            self.calc_RBtor()

        # now generate the appropriate dictionary for
        # saving the files

        tlbl = 'Time'.ljust(20) + 'Seconds   '
        rholbl = 'rho_tor'+' '*23
        lbl = 'R*Bt'+' '*20+'T*m   '
        uf_d = {'pre': 'F', 'ext': 'EQ', 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl, 'arr': self.time_u},
                         'Y': {'lbl': rholbl, 'arr': self.rho}},
                'data': {'lbl': lbl, 'arr': self.RBT}}
        ufiles.WU(uf_d, udir=self.path) 


    def plot_input(self, *argv):
        """
        Function that plots 1D (q, R.Bt) and 2D (magnetic surfaces)
        profiles. Useful for benchmark of input quantities

        """
        
        if __name__ == '__main__':
            self.pltframe = tk.TK()
        else:
            self.pltframe = tk.Toplevel()

        self.pltframe.title(str(argv[-1]))

        
        n_sp=self.nt
        if n_sp==1:
            r_nsp=1
            c_nsp=1
        elif n_sp%2==0:
            r_nsp=2
            c_nsp=n_sp/2
        else:
            r_nsp=2
            c_nsp=(n_sp+1)/2
        #fig=plt.figure()
        fig = Figure(figsize=(5,5), dpi=100)
        
        for i in range(self.nt):
            ax = fig.add_subplot(r_nsp,c_nsp,i+1)
            ax.set_title("time="+str(self.time_u[i]))
            if len(argv)==3:
                if str(argv[-1])=='OUT' or str(argv[-1]) == 'RZ LCFS after moments':
                    ax.plot(argv[0][i,:], argv[1][i,:])
                else:
                    for j in range(self.n_rho):
                        ax.plot(argv[0][i,j,:], argv[1][i,j,:])
            elif len(argv)==2:
                if str(argv[-1])=='Fpol':
                    ax.plot(np.linspace(0, 1, num=51),argv[0][i,:])
                else:
                    ax.plot(self.rho, argv[0][i,:])

        #fig=plt.gcf()
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, self.pltframe)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        button = sty.myButton(self.pltframe, 'Dismiss', self.pltframe.destroy, bgc=tksty.qtcol)
        button.pack()
        toolbar = NavigationToolbar2TkAgg(canvas, self.pltframe)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)    
            
class geo_indict(GEO):
    """
    ====
    geo_indict class (inherited from GEO class)
    ====

    This class uses the input dictionary for reading from LIUQE
    Input:
    ---------------

    A dictionary containing all the needed information. The dictionary may
    be created through a GUI or directly given
    Dictionary is of the form
    {'shot': shot, 'tbeg': tbeg, 'tend':tend,'dt':time resolution, 
    'Nmom':number of fourier moments,'Mom type':4, 
    'n_rho': number of point in toroidal rho,'n_the': number of point in theta, 
    'run': run of the shot}
    
    e.g.
    inDict={'shot':53778, 'dt':0.02, 'tbeg':1.4, 'tend':1.4, 'n_rho':41, 'n_the':101, 'n_rho_eq':40, 'path': YOUR HOME, 'run':'A01'}

    The path where the ufiles are stored is automatically defined
    $HOME/tr_client/TCV/#SHOT/#RUN
    
    Methods:
    ---------------
    _init_indict: gets the equilibrium, LCFS, defines grids
    _mag_axis: calculates position of axis
    _calc_LCFS: calculates the LCFS position
    _eqflux: calculates magnetic surfaces position
    _calc_RBtor: calculates and returns poloidal flux on appropriate grid
    _calc_qprof: calculates and returns safety factor profile on appropriate grid

    Attributes:
    ---------------
        rz_grid: r, z grid used in the computation of the field
        rho: rho gridspec
        theta: theta gridspec
    """
    def __init__(self, indict):
        """
        Initialise class using dictionary 
        """
        GEO.__init__(self, indict)
        self._init_indict()
        
    def _init_indict(self):
        print("DOING GEO INDICT")
        self.dict_flag=1
        # we have not implemented a function
        # with a variable number of moments
        #self.mom_order = self.geodict['Nmom']

        #self.path = os.path.expanduser("~") + \
        #            '/tr_client/TCV/' + str(self.shot)
        # now open the connection and define the equilibrium
        # open the appropriate equilibrium
        self.eq = eqtools.TCVLIUQETree(self.shot)
        # we have a series of things which are common
        # as the time basis and the appropriate grid
        t = self.eq.getTimeBase()
        _ind = ((t >= self.tbeg) & (t < self.tend))
        if len(np.where(_ind==True))<=3:
            _ind = ((t >= self.tbeg-0.05) & (t < self.tend+0.05))
            
        if self.tbeg==self.tend:
            t_req = np.array([self.tbeg])
        else:
            t_req = np.arange(self.tbeg, self.tend + self.dt, self.dt)

        self.time_u = t_req

        self._calc_LCFS()
        
        # this is the grid used throughout the class
        self.R = np.linspace(np.min(self.rc), np.max(self.rc), num=2048)
        self.Z = np.linspace(np.min(self.zc), np.max(self.zc), num=2048)

        # this is the time basis and the number of points
        # in time through the class
        self.nt = self.time_u.size
        self.rz_grid = self.R, self.Z
        # this is the rho and theta definition throught the class
        self.rho = np.linspace(0, 1, num=self.n_rho)
        self.theta = np.linspace(0, 2*np.pi, self.n_the, endpoint=False)


        #==========================================
        # initialising COCOS conversion factors
        self._coord_change()
        
        #self.plot_input(self.rc, self.zc, 'OUT')
        
        
    def _mag_axis(self):
        """
        Hidden method to retrieve the position of magnetic axis and interpolate
        it in the requested time samples
        """
        t = self.eq.getTimeBase()
        _ind = ((t >= self.tbeg-5*self.dt) & (t < self.tend+5*self.dt))

        _tmp_r_axis = self.eq.getMagR()[_ind]
        param_raxis = interpolate.interp1d(t[_ind],  _tmp_r_axis)
        self.r_axis = param_raxis(self.time_u)

        _tmp_z_axis = self.eq.getMagZ()[_ind]
        param_zaxis = interpolate.interp1d(t[_ind],  _tmp_z_axis)
        self.z_axis = param_zaxis(self.time_u)    


    def _calc_LCFS(self):
        """
        Function to compute the LCFS and thus the RZ limits
        read the information on equilibrium from LIUQE. We use directly
        the saved r and z contour of the LCFS already limited
        
        """
        t = self.eq.getTimeBase()
        _ind = ((t >= self.tbeg) & (t < self.tend))
        if len(np.where(_ind==True))<=4:
            _ind = ((t >= self.tbeg-0.1) & (t < self.tend+0.1))
            
        _rc = self.eq.getRLCFS()[_ind, :]
        _rc = _rc[:,~np.isnan(_rc).any(0)] #REMOVES NANs
        num_LCFS = _rc.shape[1]
        x_numLCFS = range(num_LCFS)
        #if self.tbeg!=self.tend:
        param_rc = interpolate.RectBivariateSpline(t[_ind], x_numLCFS , _rc)
        self.param_rc = param_rc
        self.rc = param_rc(self.time_u, x_numLCFS)
        #else:
        #    self.rc = _rc

        _zc = self.eq.getZLCFS()[_ind, :]
        _zc = _zc[:,~np.isnan(_zc).any(0)] #REMOVES NANs
        #if self.tbeg!=self.tend:
        param_zc = interpolate.RectBivariateSpline(t[_ind], x_numLCFS , _zc)
        self.param_zc = param_zc
        self.zc = param_zc(self.time_u, x_numLCFS)
        #else:
        #    self.zc=_zc
        #=========================================
        if self.tbeg == self.tend:
            self.rc = self.rc[0,:]
            self.zc = self.zc[0,:]


        
#    def _calc_psi_deriv(self):
#        """
#        Compute the derivative of the poloidal flux on a refined grid which
#        will be used then for computing of the radial and vertical component of
#        the magnetic field.
#        It can be done by computing on a finer grid (128x128)
#        within the separatrix
#        for each time
 #       """
 #       self.dpsidR = np.zeros((self.nt, self.R.size, self.Z.size))
 #       self.dpsidZ = np.zeros((self.nt, self.R.size, self.Z.size))
 #       for time, itime in zip(self.time_u, range(self.nt)):
 #           #psiN = self.eq.rz2psi(self.R, self.Z, time, make_grid=True)
 #           psiN = self.eq.rz2psi(self.R, self.Z, time, make_grid=True)
 #           #psiN = self.coord_change_psi(psiN_t)
#            deriv = np.gradient(psiN)
#            # Note np.gradient gives y
#            # derivative first, then x derivative
 #           ddR = deriv[1]
#            # ddR = self.psi(Rgrid,Zgrid,dx=1)
#            ddZ = deriv[0]
#            # ddZ = self.psi(Rgrid,Zgrid,dy=1)
#            dRdi = 1.0/np.gradient(self.R)
#            dZdi = 1.0/np.gradient(self.Z)
#            self.dpsidR[itime, :, :] = ddR*dRdi[np.newaxis, :]
#            self.dpsidZ[itime, :, :] = ddZ*dZdi[:, np.newaxis]

#    def calc_bfield(self):
#        """
#        Calculate the radial and vertical component of the
#        poloidal magnetic field in a prescribed grid
#        These are the variables BRRZ and BZRZ of the plasma_state
#
#        """
#        self._calc_psi_deriv()
#        self.BRRZ = -1.0*self.dpsidZ/self.R[np.newaxis, :]
#        self.BZRZ = self.dpsidR/self.R[np.newaxis, :]
#        return self.BRRZ, self.BZRZ

        
#    def calc_area(self):
#        """
#        Calculate the area of the surface embedded in toroidal
#        iso-flux surfaces. It use the approximate r_plot, z_plot
#        from the computation of eqflux so that we are sure the points are
#        continous
#
#        """
#        try:
#            self.moments.mean()
#        except:
#            self._eqflux()
#
#        self.area = np.zeros((self.nt, self.n_rho))
#        for itime in range(self.nt):
#            for jrho in range(self.n_rho):
#                _dummy = greenArea(self.r_plot[itime, jrho, :],
#                                   self.z_plot[itime, jrho, :])
#                self.area[itime, jrho] = np.abs(_dummy)
#        return self.area

    def _eqflux(self):
        """
        It compute the appropriate toroidal flux in R, Z and the appropriate
        momentum flux representation. It uses eqtools
        to move from psiRZ to psitorRZ
        and then loop for the appropriate momenta.
        In this way we have appropriate values
        """
        #check if the magnetic axis position has been computed, otherwise compute it
        try:
            self.r_axis.mean()
        except:
            self._mag_axis()

        # ok now we can compute the appropriate momentum
        # we need on a rhotor grid
        self.moments = np.zeros((self.nt, self.n_rho,
                                 self.mom_order, self.mom_type))
        # the moments are a function of
        # (time, surface label, mom_order, mom_type)
        # R (time, rho, theta) this will be saved in ufile
        self.r_plot = np.zeros((self.nt, self.n_rho, self.n_the))
        # Z (time, rho, theta) this will be saved in ufile
        self.z_plot = np.zeros((self.nt, self.n_rho, self.n_the))
        self.cont = np.zeros((self.nt))
        self.tempz = np.zeros((self.nt,self.n_rho))
        # now we have to iterate on the time and
        for tshot, itime in zip(self.time_u, range(self.nt)):
            print('TSHOT', tshot)
            # for each time the grid is chosen within the LCFS
            torFlux_t = self.eq.rz2phinorm(self.R, self.Z, tshot, make_grid=True, sqrt=True)
            #torFlux = self.coord_change_torflux(torFlux_t)
            torFlux = torFlux_t
            # now we can build the appropriate contour
            rcont = contourc.contourc(self.R, self.Z, torFlux, self.n_rho)
            # we limit to the point in rho<1 afterwards we will add the
            # computed separatrix moment representation
            for jrho in np.linspace(1,
                                    self.n_rho-2,
                                    num=self.n_rho-2, dtype='int'):
                _rrzz = np.asarray(rcont[jrho])
                if _rrzz.ndim > 1:
                    rrzz = np.squeeze(_rrzz)
                else:
                    rrzz = np.vstack([_rrzz[i] for i in range(_rrzz.size)])
                
                arr = descu(rrzz[:, 0], rrzz[:, 1], self.mom_order)
                self.moments[itime, jrho, :, 0] = arr[:, 0]
                self.moments[itime, jrho, :, 1] = arr[:, 1]
                self.moments[itime, jrho, :, 2] = arr[:, 2]
                self.moments[itime, jrho, :, 3] = arr[:, 3]
                self.r_plot[itime, jrho, :], self.z_plot[itime, jrho, :] = mom2rz.mom2rz(arr[:, 0], arr[:, 1],
                                                                                         arr[:, 2], arr[:, 3], 
                                                                                         nthe=self.n_the,
                                                                                         endpoint = True)
	        # now add the same evaluation for the separatrix assumed at rho=1
            if self.tbeg!=self.tend:
                xd, yd = self.rc[itime,: ], self.zc[itime,: ]
            else:
                xd, yd = self.rc, self.zc
            arr = descu(xd[~np.isnan(xd)],yd[~np.isnan(xd)], self.mom_order)
            self.moments[itime, -1, :, 0] = arr[:, 0]
            self.moments[itime, -1, :, 1] = arr[:, 1]
            self.moments[itime, -1, :, 2] = arr[:, 2]
            self.moments[itime, -1, :, 3] = arr[:, 3]
            self.r_plot[itime, -1, :], self.z_plot[itime, -1, :] = mom2rz.mom2rz(arr[:, 0], arr[:, 1],
                                                                                 arr[:, 2], arr[:, 3], 
                                                                                 nthe=self.n_the,
                                                                                 endpoint = True)
            # now add the magnetic axis (rho=0)
            r0=np.full(self.n_the, self.r_axis[itime])
            z0=np.full(self.n_the, self.z_axis[itime])
            self.r_plot[itime, 0, :] = r0
            self.z_plot[itime, 0, :] = z0

	
	#self.plot_input(self.r_plot[:,-1,:], self.z_plot[:,-1,:], 'RZ LCFS after moments')

               
    def _calc_RBtor(self):
        """
        Hidden method to compute the Fprofile on the appropriate rho-toroidal
        grid 

        """
        # load the q values as saved in the MDSplus Tree
        t = self.eq.getTimeBase()
        # define a fake rhop
        rhop_value = np.linspace(0, 1, num=51, dtype='float')
        
        _indx = ((t >= self.tbeg-5*self.dt) & (t < self.tend+5*self.dt))

        # the profile is defined in an equidistant poloidal grid
        # which we need to convert in the appropriate toroidal grid
        _temp_fprof = self.eq.getF()[_indx,:]
        #interpolation of f_prof in the right time values
        param_f_time = interpolate.RectBivariateSpline(t[_indx], np.linspace(0,1,_temp_fprof.shape[1]), _temp_fprof)
        self._fprof = param_f_time(self.time_u, rhop_value)

        #conversion of rho_pol in rho_tor(rho_pol,time)
        _rho_tor_temp = self.eq.psinorm2phinorm(rhop_value, t[_indx], each_t=True, sqrt=True)
        #interpolation of rho_tor on desired time
        param_rhotor_time = interpolate.RectBivariateSpline(t[_indx], rhop_value, _rho_tor_temp)
        self._rho_tor = param_rhotor_time(self.time_u, rhop_value)
        

        # now we have the f profile defined on self.time_u and _rho_tor, so we interpolate in order to
        # have it on self.time_u and self.rho
        self.RBT = np.ones(((self.nt), (self.n_rho)))
        #slice_rho_tor = np.zeros(self.n_rho)
        slice_f = np.zeros(self.n_rho)
        for i in range(self.nt):
            slice_rhotor = self._rho_tor[i,:]
            slice_f = self._fprof[i,:]
            #slice_f = slice_f[::-1]
            param_f_rho = interpolate.interp1d(slice_rhotor, slice_f)
            self.RBT[i,1:-1] = param_f_rho(self.rho[1:-1])

            #now add last point, which couldn't be included before:
            lastpoint = slice_f[-1]
            self.RBT[i,-1]=lastpoint

            #now add first point
            firstpoint = slice_f[0]
            self.RBT[i,0]=firstpoint

        #self.RBT = self.RBT*1.24

        #self.RBT = np.linspace(1.25, 1.26, self.n_rho)
        #self.RBT = self.RBT[np.newaxis,:]
        #print('SORTING')
        #for i in range(self.nt):
        #    self.RBT[i,:] = np.sort(self.RBT[i,:])

        #self.coord_change_F()

        self.RBT = np.absolute(self.RBT)
        

        

    def _calc_qprof(self):
        """
        Hidden method to compute the qprofile on the appropriate rho-toroidal
        grid 
        the q got from eqtools is defined on PSIGRID

        """
        # load the q values as saved in the MDSplus Tree
        t = self.eq.getTimeBase()
        # define a fake psigrid
        psiFake = np.linspace(0, 1, num=51)
        _indx = ((t >= self.tbeg-5*self.dt) & (t < self.tend+5*self.dt))

        # the profile is defined POLOIDAL FLUX GRID
        # which we need to convert in the appropriate toroidal grid
        _temp_qprof = self.eq.getQProfile()[_indx, :]
        #interpolation of q_prof in the right time values
        param_q_time = interpolate.RectBivariateSpline(t[_indx], np.linspace(0,1,_temp_qprof.shape[1]), _temp_qprof)
        _qprof = param_q_time(self.time_u, psiFake)

        #conversion of rho_pol in rho_tor
        psiFake = np.linspace(0, 1, num=51, dtype='float')
        #psinorm2phinorm gets the FLUX, not rho (from docs), but if sqrt=True returns rho
        _rho_tor_temp = self.eq.psinorm2phinorm(psiFake, t[_indx], each_t=True, sqrt=True)
        #interpolation of rho_tor on desired time
        param_rhotor_time = interpolate.RectBivariateSpline(t[_indx], psiFake, _rho_tor_temp)
        _rho_tor = param_rhotor_time(self.time_u, psiFake)

#        plt.plot(_rho_tor[0,:], _qprof[0,:])
        
        # now we have the q profile defined on self.time_u and _rho_tor, so we interpolate in order to
        # have it on self.time_u and self.rho
        self.qprof = np.zeros(((self.nt), (self.n_rho)))
        slice_q = np.zeros(np.size(psiFake))
        for i in range(self.nt):
            slice_rhotor = _rho_tor[i,:]
            slice_q = _qprof[i,:]
            param_q_rho = interpolate.interp1d(slice_rhotor, slice_q)
            self.qprof[i,1:-1] = param_q_rho(self.rho[1:-1])

            #now add last point, which couldn't be included before:
            lastpoint = slice_q[-1]
            self.qprof[i,-1]=lastpoint

            #now add first point
            firstpoint = slice_q[1]
            self.qprof[i,0]=firstpoint

         
        #self.coord_change_q()
        #self.plot_input(self._qprof, 'Q')
            


class geo_eqdsk(GEO):
    """
    ====
    geo_eqdsk class (inherited from GEO class)
    ====

    This class uses the eqdsk file for storing equilibrium
    Input:
    ---------------
    (1) A dictionary containing all the needed information. The dictionary may
    be created through a GUI or directly given
    Dictionary is of the form
    {'shot': shot, 'tbeg': tbeg, 'tend':tend,'dt':time resolution, 
    'Nmom':number of fourier moments,'Mom type':4, 
    'n_rho': number of point in toroidal rho,'n_the': number of point in theta, 
    'run': run of the shot}
    e.g.
    inDict={'shot':53778, 'dt':0.02, 'tbeg':1.4, 'tend':1.4, 'n_rho':41, 'n_the':101, 'n_rho_eq':40, 'path': YOUR HOME, 'run':'A01'}

    (2) eqdsk file to use 

    
    Methods:
    ---------------
    _init_eqdsk: gets the data from eqdsk file (axis, LCFS, flux in grid and profile)
    _psigrid2phigrid: converts from poloidal flux to toroidal flux (needed for e.g. nubeam)
    _eqflux: calculates magnetic surfaces position


    Attributes:
    ---------------
        rz_grid: r, z grid used in the computation of the field
        rho: rho gridspec
        theta: theta gridspec
    """

    
    def __init__(self, indict,infile_eqdsk):
        GEO.__init__(self, indict)
        self._init_eqdsk(infile_eqdsk)
        
    def _init_eqdsk(self, infile_eqdsk):
        """
        function for import from eqdsk file
        this is the structure of the eqdsk struct:
        
            self.comment=comment
            self.switch=switch
            self.nrbox=nrbox
            self.nzbox=nzbox
            self.rboxlength=rboxlength
            self.zboxlength=zboxlength
            self.R0EXP=R0EXP
            self.rboxleft=rboxleft
            self.Raxis=Raxis
            self.Zaxis=Zaxis
            self.psiaxis=psiaxis
            self.psiedge=psiedge
            self.B0EXP=B0EXP
            self.Ip=Ip
            self.T=T
            self.p=p
            self.TTprime=TTprime
            self.pprime=pprime
            self.psi=psi
            self.q=q
            self.nLCFS=nLCFS
            self.nlimits=nlimits
            self.R=R
            self.Z=Z
            self.R_limits=R_limits
            self.Z_limits=Z_limits
            self.R_grid=R_grid
            self.Z_grid=Z_grid
            self.psi_grid=psi_grid  ### POLOIDAL FLUX
            self.rhopsi=rhopsi
        
        """
        self.theta = np.linspace(0, 2*np.pi, self.n_the, endpoint=False)
        self.eqdsk_flag=1
        self.time_u = np.array([self.tbeg])
        self.eqdsk= ReadEQDSK_MV.ReadEQDSK(infile_eqdsk)
        self.nrho_eqdsk = len(self.eqdsk.psi_grid)
        self.psi_grid = (self.eqdsk.psi-self.eqdsk.psiaxis)/(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        print (self.psi_grid.shape)
        rho_eqdsk = self.eqdsk.rhopsi
        # Conversion of rho in the right grid. rhopsi is the sqrt of the flux
        #self.rho = rho_eqdsk
        self.fakerho = np.linspace(0,1,self.n_rho)
        # poloidal normalised
        psi_t = self.eqdsk.psi_grid
        psi_t = (psi_t-self.eqdsk.psiaxis)/(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        self.param_psi = interpolate.interp1d(rho_eqdsk, psi_t)

        #FIND LCFS
        #psi_grid = (self.eqdsk.psi-self.eqdsk.psiaxis)/(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        #cs=plt.contour(self.eqdsk.R_grid, self.eqdsk.Z_grid, psi_grid, 1)
        #tmp=cs.collections[0].get_paths()[0]
        self.rc, self.zc = np.array(self.eqdsk.R), np.array(self.eqdsk.Z)

        qprof_t = self.eqdsk.q
        self.param_q = interpolate.interp1d(rho_eqdsk, qprof_t)
        self._psi2phi()
        
        # poloidal flux (F in grad-shafranov equation)
        RBT_t = self.eqdsk.T 
        self.param_RBT = interpolate.interp1d(rho_eqdsk, RBT_t)
        RBT = (-1)*self.param_RBT(self.rho)
        self.RBT = np.expand_dims(RBT, axis=0)                      

        # q profile
        qprof = self.param_q(self.rho)
        self.qprof = np.expand_dims(qprof, axis=0)

        # magn. axis
        self.r_axis = self.eqdsk.Raxis
        self.z_axis = self.eqdsk.Zaxis

        #R and Z for spacing
        self.R = np.linspace(np.min(self.rc), np.max(self.rc), num=4096)
        self.Z = np.linspace(np.min(self.zc), np.max(self.zc), num=4096)

        #defining number of timesteps (for _eqflux)
        self.nt = len(self.time_u)
        
    def _psi2phi(self):
        """
        Converts psi 2 phi
        """
        tmpnum=1000
        #self.phi = np.zeros(tmpnum)
        locq   = self.param_q(np.linspace(0,1,tmpnum)) #augmenting precision near the core
#        locpsi = self.param_psi(np.linspace(0,1,tmpnum))*(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        locpsi = np.linspace(0,1,tmpnum)*(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        phi = integrate.cumtrapz(locq,locpsi)
        phi = np.concatenate([[0], phi])
        self.phi_edge = np.max(np.abs(phi))
        #phi = phi/np.max(phi)
        self.param_phi = interpolate.interp1d(np.linspace(0,1,tmpnum), phi) #interpolate on psinorm grid
        self.phi = self.param_phi(self.fakerho) # compute on fake rho grid
        self.phinorm = self.phi/self.phi_edge
        self.rhotor = np.abs(self.phinorm)**0.5 #defined on fake rho grid
        self.rho = self.rhotor
        #plt.plot(locpsi, locq, label='q')
        #plt.plot(locpsi, phi, label='phi')

    def _psigrid2phigrid(self):
        """
        Converts between the two grids just by using phi(psi)
        """
        try:
            self.phi_edge.mean()
        except:
            self._psi2phi()
        #interpolate on psi
        print(self.eqdsk.psi_grid.shape)
        param_psi = interpolate.interp2d(self.eqdsk.R_grid,self.eqdsk.Z_grid, self.psi_grid)
        self.psi_grid = param_psi(self.R, self.Z)
        
        tmpnum=1000
        locpsi = np.linspace(0,1,tmpnum)
        locphi = self.param_phi(locpsi)/self.phi_edge #more fine phi grid on psi grid
        #plt.plot(locpsi, locphi,'k')
        #plt.plot(self.psi, self.phinorm, 'r')
        phi_grid = np.zeros(np.shape(self.psi_grid))
        for index,psi_t in np.ndenumerate(self.psi_grid):
            ind = np.abs(locpsi-psi_t).argmin()
            phi_grid[index] = locphi[ind]
#        phi_grid=phi_grid/np.max(phi_grid)
        self.phi_gridnosmoot = phi_grid
        #param_phi=interpolate.interp2d(self.eqdsk.R_grid,self.eqdsk.Z_grid, phi_grid)
        #self.phi_grid = param_phi(self.R, self.Z)
        self.phi_grid = self.phi_gridnosmoot
        self.rho_tor_grid = np.sqrt(np.abs(self.phi_grid))
                    
    def _eqflux(self):
        """
        THIS IS FOR EQDSK FILES
        It compute the appropriate toroidal flux in R, Z and the appropriate
        momentum flux representation. 
        it uses _psigrid2phigrid for moving btw psi and phi
        and then loop for the appropriate momenta.
        """
        try:
            self.phi_grid.mean()
        except:
            self._psigrid2phigrid()
        # ok now we can compute the appropriate momentum
        # we need on a rhotor grid
        self.moments = np.zeros((self.nt, self.n_rho,
                                 self.mom_order, self.mom_type))
        # the moments are a function of
        # (time, surface label, mom_order, mom_type)
        # R (time, rho, theta) this will be saved in ufile
        self.r_plot = np.zeros((self.nt, self.n_rho, self.n_the))
        # Z (time, rho, theta) this will be saved in ufile
        self.z_plot = np.zeros((self.nt, self.n_rho, self.n_the))
        self.cont = np.zeros((self.nt))
        self.tempz = np.zeros((self.nt,self.n_rho))

        torFlux = self.rho_tor_grid
        # now we can build the appropriate contour
        rcont = contourc.contourc(self.R, self.Z, torFlux, self.n_rho)
        # we limit to the point in rho<1 afterwards we will add the
        # computed separatrix moment representation
        for jrho in np.linspace(1,
                                self.n_rho-2,
                                num=self.n_rho-2, dtype='int'):
            _rrzz = np.asarray(rcont[jrho])
            if _rrzz.ndim > 1:
                rrzz = np.squeeze(_rrzz)
            else:
                rrzz = np.vstack([_rrzz[i] for i in range(_rrzz.size)])
                
            arr = descu(rrzz[:, 0], rrzz[:, 1], self.mom_order)
            self.moments[0, jrho, :, 0] = arr[:, 0]
            self.moments[0, jrho, :, 1] = arr[:, 1]
            self.moments[0, jrho, :, 2] = arr[:, 2]
            self.moments[0, jrho, :, 3] = arr[:, 3]
            self.r_plot[0, jrho, :], self.z_plot[0, jrho, :] = mom2rz.mom2rz(arr[:, 0], arr[:, 1],
                                                                             arr[:, 2], arr[:, 3], 
                                                                             nthe=self.n_the,
                                                                             endpoint = True)

	# now add the same evaluation for the separatrix assumed at rho=1
        xd, yd = self.rc, self.zc
        arr = descu(xd[~np.isnan(xd)],yd[~np.isnan(xd)], self.mom_order)
        self.moments[0, -1, :, 0] = arr[:, 0]
        self.moments[0, -1, :, 1] = arr[:, 1]
        self.moments[0, -1, :, 2] = arr[:, 2]
        self.moments[0, -1, :, 3] = arr[:, 3]
        self.r_plot[0, -1, :], self.z_plot[0, -1, :] = mom2rz.mom2rz(arr[:, 0], arr[:, 1],
                                                                     arr[:, 2], arr[:, 3], 
                                                                     nthe=self.n_the,
                                                                     endpoint = True)
        # now add the magnetic axis (rho=0)
        r0=np.full(self.n_the, self.r_axis)
        z0=np.full(self.n_the, self.z_axis)
        self.r_plot[0, 0, :] = r0
        self.z_plot[0, 0, :] = z0

	
	#self.plot_input(self.r_plot[:,-1,:], self.z_plot[:,-1,:], 'RZ LCFS after moments')



