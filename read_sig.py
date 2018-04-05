from __future__ import print_function
import os
import ufiles
import tcv
import MDSplus as mds
import eqtools
import numpy as np
from scipy import interpolate
import Tkinter as tk
import tksty

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

tlbl = 'Time'.ljust(20) + 'Seconds   '
ylbl = 'rho_tor' + ' ' * 23

col = 2 * ['#c00000', '#00c000', '#0000c0', '#b0b000', '#b000b0', '#00b0b0']
fsize = 14
sty = tksty.TKSTY()


# this is the class for reading the 1D signals and eventually store the
# corresponding ufiles
class DWR1:
    """
    Class for reading the 1D signals needed for the TRANSP simulation
    It uses the same dictionary as the GEO but we limit to the
    shot number,tbeg and tend. The signals definition is embedded
    into the class

    dictIn= {'shot': shot, 'tbeg': tbeg, 'tend':tend,
    'dt': time resolution, 'path': path where to store the outputs}

    ---------------
    Methods:
    store_1d = Store the 1D of the signals (Ip,BTF,ZEFF,VLOOP)
    read_1d = gives back a dictionary of the type
    {'IP':{'time':time,'data':data}} with the signals
    ---------------
    Attributes:

    """

    def __init__(self, indict):
        
        self.shot = indict['shot']
        self.tbeg = indict['tbeg']
        self.tend = indict['tend']
        self.dt = indict['dt']
        # open the connection to tcvdata
        self.conn = tcv.shot(self.shot)
        self.path = indict['path']
        self.tree = mds.Tree('tcv_shot', self.shot)
        #self.path = os.path.expanduser("~") + '/tr_client/TCV/' + str(self.shot)
        # now check existence of path otherwise creates it
        if os.path.exists(self.path) is False:
            os.system('mkdir -p %s' % self.path)

        self._univec = {}  # this is a dictionary of object
        self.time_u = np.arange(self.tbeg, self.tend, self.dt)

        # this is the dictionary used for defining the proper
        # 1D signals which will be read
        
#        TRANSP DATA NEEDED (1D)
#        data=struct('IPL',{},'RBTOR',{},'NE',{},'TE',{},'TI',{},'PRAD',{}, 'Q',{}, ...
#        'VTOR',{},'LI',{},'WMHD',{},'VLOOP',{},'BETAT',{},'NEUT',{},'OH',{},'SXR',{});
        
        self.signals = {'IP': {'suff': 'CUR', 'prefix': 'I',
                               'string': r"\results::I_P",
                               'yl': r'I$_{pl}$',
                               'ytk': [0, 1e5, 2e5, 3e5, 4e5],
                               'lbl': 'Plasma Current Amps'},
                        'BTF': {'suff': 'RBT', 'prefix': 'B',
                                'string':r"\magnetics::RBPHI",
                                'yl': 'R*B$_{o}$',
                                'ytk': [300., 350, 400., 450],
                                'lbl': 'Rp*Bt               T.cm'},
                        'ZEFF': {'string': r'\results::conf:z_eff',
                                 'yl': r'Z$_{eff}$',
                                 'ytk': [1, 2, 3], 'lbl': 'Zeff  ',
                                 'suff': 'EFF', 'prefix': 'Z'},
                        'VLOOP': {'string': r"\magnetics::vloop",
                                  'yl': r'V loop [V]',
                                  'ytk': [-2, 0, 2, 4],
                                  'lbl': 'Vloop'.ljust(20) + '[V]',
                                  'suff': 'LOOP', 'prefix': 'V'}}



        print("\n")
        print("===================")
        print("Initialisation of 1D signals Done")
        print("===================")
        print("\n")



    def _readJB(self):
        """
        This is needed in order to write to the namelist the values
        in a fast way
        """
        try:
            return np.sign(np.mean(self.rsig['IP']['data'])),\
        np.sig(np.mean(self.rsig['BTF']['data']))
        except:
            self.signs={'IP':0, 'BTF':0}

            stri = self.signals['IP']['string']
            data = self.conn.tdi(stri).values
            self.signs['IP']=np.sign(np.mean(data))
                
            stri = r'\tcv_shot::top.magnetics.measurements.currents:rbphi'
            data = self.tree.getNode(stri).data()
            self.signs['BTF']=np.sign(np.mean(data))

            return int(self.signs['IP']), int(self.signs['BTF'])

    def writeJB(self):
        """
        Writes file with the signs of J and B           
        """
        try:
            self.signs['IP']
        except:
            self._read_JB()
            
        fname=self.path+'BTIP_SIGN.DAT'
        outfile = open(fname, 'w')
        outfile.write('BT sign:'+str(self.signs['BTF'])+'\n')
        outfile.write('IP sign:'+str(self.signs['IP']))
        
        outfile.close()
        
    def read_1d(self):
        """
        Returns the dictionary of the signals (Ip,RBT,ZEFF,VLOPP,NBI)
        """
        try:
            self.rsig
        except:
            self._read_1d()
        print("\n")
        print("===================")
        print("END READING 1D")
        print("===================")
        print("\n")
       
        
    def _read_1d(self):

        """
        We define the appropriate attribute rsig, which is a dictionary
        with the signals and the time at the resolution needed
        """
        # check we define the appropriate
        # dictionary of UnivariateSpline
        # representation
        try:
            self._univec()
        except:
            self._getUnivecSpline()
            
        self.rsig = {}
        for k in self._univec.keys():
            self.rsig[k] = dict([('data',
                                  self._univec[k]['spline'](self.time_u)),
                                 ('time', self.time_u)])
            
    def _getUnivecSpline(self):

        """
        Return a dictionary of UnivecSpline representation
        of the signals so that any time resolution can be obtained
        """

        # define the new time basis which will be the real
        # time basis of the signals saved
        if self.tbeg==self.tend:
            self.time_u = np.array([self.tbeg])
            
        for k in self.signals.keys():
            stri = self.signals[k]['string']
            print('Reading signal '+stri)
            # now read the signals
            try:
                data = self.conn.tdi(stri).values
                tim = data.dim_0.values
            except:
                data = self.tree.getNode(stri).data()
                if k == 'BTF':
                    data = data[0,:]
                    data = data*1e2 #Tm to Tcm 
                if k == 'VLOOP':
                    data = data[0,:]
                tim = self.tree.getNode(stri).getDimensionAt(0).data()
            
            # limit to the times between tbeg and tend
            _iidx = ((tim >= self.tbeg-self.dt) & (tim <= self.tend+self.dt))
            if len(tim[_iidx])<=4:
                _iidx = ((tim >= self.tbeg-0.1) & (tim <= self.tend+0.1))
            dat = data[_iidx]
            if k == 'IP':
                dat = np.fabs(dat)
            tim = tim[_iidx]
            dummy = interpolate.InterpolatedUnivariateSpline(tim, \
                                                             dat, ext=0)
            self._univec[k] = dict([('spline', dummy)])

        print("\n")
        print("===================")
        print("Spline of 1D data done")
        print("===================")
        print("\n")

    def set_zeff(self):
        """
        Allows to set a constant value of Zeff. Needed because I don'trust
        the values by diagnostic
        """
        try:
            self.rsig['ZEFF']
        except:
            self._read_1d()
            
        zeff = input(' Insert value of Zeff to use (constant over time) \n')
        self.zeff=np.full(len(self.time_u), zeff, dtype=float)
        x,y = self.time_u, self.zeff
        dummy = interpolate.InterpolatedUnivariateSpline(x, y, \
                                                         ext=0)
        
        self._univec['ZEFF'] = dict([('spline', dummy)])
        self.rsig['ZEFF'] = dict([('data',
                                  self._univec['ZEFF']['spline'](self.time_u)),
                                 ('time', self.time_u)])
        
        
    def store_1d(self):
        try:
            self.rsig
        except:
            self._read_1d()

        self._store_1d()        
        
    def _store_1d(self):
        """
        This method store the data which has already been read. It does not
        store the ECRH which have multiple sources
        """

        # first of all check the existence of self.rsig
        try:
            self.rsig
        except:
            self._read_1d()

        # now we cycle to save the appropriate ufiles
        for k in self.signals.keys():
            uf_d = {'pre': self.signals[k]['prefix'],
                    'ext': self.signals[k]['suff'],
                    'shot': self.shot,
                    'grid': {'X': {'lbl': tlbl,
                                   'arr': self.rsig[k]['time']}},
                    'data': {'lbl': self.signals[k]['lbl'],
                             'arr': self.rsig[k]['data']}}
            
            ufiles.WU(uf_d, udir=self.path)



    def plot_input(self):
        """
        Function that plots 1D quantities. 
        Useful for benchmark of input quantities

        """
#        if __name__ == '__main__':
#            self.pltframe = tk.TK()
#        else:
#            self.pltframe = tk.Toplevel()
#
#        self.pltframe.title('1D')

        
        #n_sp=self.nt
        #if n_sp%2==0:
        #    r_nsp=2
        #    c_nsp=n_sp/2
        #else:
        #    r_nsp=2
        #    c_nsp=(n_sp+1)/2

        r_nsp=2; c_nsp=3
        #fig=plt.figure()
        fig = plt.figure(figsize=(15,10))
        key_arr =  [ 'IP', 'BTF', 'ZEFF', 'VLOOP']
        key_arr += [ 'NBI', 'EC']

        for indi, i in enumerate(key_arr):
            ax = fig.add_subplot(r_nsp, c_nsp,indi+1)
            try:
                name = self.signals.keys()[indi]
                data = self.rsig[name]['data']
                ax.plot(self.time_u, data, 'k')
            except:
                name = self.chann_sig.keys()[indi-len(self.signals.keys())]
                for jj in self.rsig_channel[name].keys():
                    data = self.rsig_channel[name][jj]['data']
                    ax.plot(self.time_u, data, label=jj)

            title = name
            ax.set_title(title)
            ax.legend(loc='best')
        fig.tight_layout() 
        plt.show()
        #fig=plt.gcf()
#        canvas = FigureCanvasTkAgg(fig, self.pltframe)
#        canvas.show()
#        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#        button = sty.myButton(self.pltframe, 'Dismiss', self.pltframe.destroy, bgc=tksty.qtcol)
#        button.pack()
        

        
        
class DWR2:
    
    def __init__(self, indict):
        
        self.indict = indict
        self.shot = indict['shot']
        self.tbeg = indict['tbeg']
        self.tend = indict['tend']
        self.dt = indict['dt']
        self.n_rho = indict['n_rho']
        self.n_rho_eq = indict['n_rho_eq']
        
        if self.tbeg==self.tend:
            self.time_u = np.array([self.tbeg])
        else:
            self.time_u = np.arange(self.tbeg, self.tend, self.dt)
            
        self.rho = np.linspace(0, 1, num=self.n_rho)
        # open the tree
        self.tree = mds.Tree('tcv_shot', self.shot)
        # open the call to eqtools for handling properly the
        # coordinate transformation
        self.eq = eqtools.TCVLIUQETree(self.shot)
        # path for saving the UFILEs
        self.path = indict['path']
        # now check existence of path otherwise creates it
        if os.path.exists(self.path) is False:
            os.system('mkdir -p %s' % self.path)
        # we build the appropriate dictionary similarly to what done
        # for the 1D signal
        self.signals = {'ne': {'string': r'\tcv_shot::top.results.thomson.profiles.auto:ne',
                               'zl': r'n$_e$ [10$^{19}$m$^{-3}$]',
                               'yl': r'$\rho_{\phi}$',
                               'suff': 'ELE', 'prefix': 'N',
                               'xlbl': tlbl, 'ylbl': ylbl,
                               'dlbl': 'NE                  [cm^-3]  '},
                        'te': {'string': r'\tcv_shot::top.results.thomson.profiles.auto:te',      
                               'zl': r'T$_e$ [eV]',
                               'yl': r'$\rho_{\phi}$',
                               'suff': 'ELE', 'prefix': 'T',
                               'xlbl': tlbl, 'ylbl': ylbl,
                               'dlbl': 'TE'.ljust(26) + '[eV]'},
                        'ti': {'string': r'\results::conf:ti',
                               'zl': r'T$_i$ [eV]',
                               'yl': r'$\rho_{\phi}$',
                               'suff': 'ION', 'prefix': 'T',
                               'xlbl': tlbl,
                               'ylbl': ylbl,
                               'dlbl': 'TI'.ljust(26) + '[eV]'}}

        print("\n")
        print("===================")
        print("Initialisation of 2D signals  Done")
        print("===================")
        print("\n")

    def _quit(self):
       self.fmframe.destroy()
       self.pltframe.destroy()
        
    def _smooth(self,x,y,w_len,beta):
        """
        Hidden method to smooth a signal (like pe)
        """
        if w_len%2 != 0:
            half_win = (w_len - 1)/2
        else:
            half_win = w_len/2
        xx = np.r_[y[w_len-1:0:-1],y,y[-1:-w_len:-1]]
	w  = np.kaiser(w_len,beta)
	yy = np.convolve(w/w.sum(),xx,mode='valid')

        return yy[half_win:-half_win]

        
    def _getBivecSpline(self):
        """
       Hidden method for reading the signal storing their bivacspline
        representation on a grid (time, rho_tor)
        """
        # Each of the signals is returned in a rho poloidal grid
        # we must convert it for each time in the corresponding rhotor
        # and then we spline on the given number of point in rho
        self._brep = {}

        for k in self.signals.keys():
            # now read the signals
            print('Reading signal ' + self.signals[k]['string'])
            tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()
            
            # choose only the timining between tbeg and tend
            if self.tbeg!=self.tend:
                if k=='vtor':
                    #increase time limit due to lower resolution
                    _idx=((tim >= self.tbeg-10*self.dt) & (tim <= self.tend+10*self.dt))
                else:
                    _idx = ((tim > self.tbeg) & (tim <= self.tend))
            else:
                _idx = np.argmin(tim-self.tbeg > 0)

            if k=='ti':
            # with thomson
                data = self.tree.getNode(self.signals[k]['string']).data()[_idx,:]
            else:
                data = self.tree.getNode(self.signals[k]['string']).data()[:, _idx]

            if k=='ne':
                data=data*1e-6

            tim = tim[_idx]

            
            # now the fake rhotoroidal equidistant grid
            _rhotor = np.linspace(1e-9, 1, num=128)
            # different approach if for the toroidal velocity which is not
            # saved on the conf tree
            if k == 'vtor':
                rhop = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()[_idx, :]
                #Define a specific _rhotor, otherwise it breaks at brep
                # convert rhop to rhotor for each time basis
                rhot = np.zeros((tim.size, rhop.shape[1]))
                for t, i in zip(tim, range(tim.size)):
                    rhot[i, :] = self.eq.psinorm2phinorm(rhop, t,
                                                         sqrt=True)
            elif k == 'ti':
                rhot = self.tree.getNode(r'\results::conf:rhotor').data()[_idx,:]
                rhot = np.sqrt(rhot)
            else:
                rhop = self.tree.getNode(r'\results::thomson.profiles.auto:rho').data()
                rhot = np.zeros((tim.size, rhop.shape[0]))
                for t, i in zip(tim, range(tim.size)):
                    rhot[i, :] = self.eq.psinorm2phinorm(rhop, t,
                                                         sqrt=True)

            signal = np.zeros((tim.size, _rhotor.size))
                
            for idx in range(tim.size):
                x = rhot[idx, :]
                if k == 'ti':
                    y = data[idx, :] #conf
                else:
                    y = data[:, idx] #thomson
                s_guess = np.var(y)*np.size(y)
                u,sm = interpolate.splprep([y], s=s_guess*0)
                y = interpolate.splev(sm,u)
                dummy = interpolate.interp1d(x, y) #, fill_value='extrapolate')
                # now we can have the bivariate
                signal[idx, :] = dummy(_rhotor)                               
            box = [tim.min(),tim.max(), 0, 1]
            brep = interpolate.RectBivariateSpline(tim, _rhotor, signal, bbox=box)

                
            # we thus need to transform in rhotor and then
            # spline on the appropriate grid
            self._brep[k] = dict([('spline', brep)])

    def _read_ticonf(self):
        """
        Function to read the ti data from conf nodes, which has different structure from the thomson
        """
        tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()
        # choose only the timining between tbeg and tend
        if self.tbeg!=self.tend:
                _idx = ((tim > self.tbeg) & (tim <= self.tend))
        else:
                _idx = np.argmin(tim-self.tbeg > 0)
            
    def read_2d(self):
        """
        Method to get the signal defined in an attribute of the
        class with the appropriate resolution in rho_toroidal
        and time. It create the attribute self.rsig with a dictionary
        with all the signals with corresponding time and rho basis
        """
        try:
            self._brep
        except:
            self._getBivecSpline()
            
        self.rsig = {}
        if self.tbeg != self.tend:
            time = np.arange(self.tbeg, self.tend+self.dt, self.dt)
        else:
            time=np.array([self.tbeg])
        rho_eq = np.linspace(0, 1, num=self.n_rho_eq)
        #THERE ARE TWO RHO, one for the magn eq and one for profiles.
        #rho_noteq is the one for the profiles
        rho_noteq = np.linspace(0, 1, num=self.n_rho)
        for k in self.signals.keys():
            if k != 'p':
                rho = rho_noteq
            else:
                rho = rho_eq

            try:
                if k!='p':
                    self.rsig[k] = dict([('signal',
                                          self._brep[k]['spline'](time, rho)),
                                         ('time', time),
                                         ('rho', rho)])
                else:
                    self.rsig[k] = dict([('signal',
                                          self._brep[k]['spline'](time, rho)),
                                         ('time', time),
                                         ('rho', rho)])                    

            except:
               self.rsig[k]=dict([('signal',0), ('time', 0),('rho', 0)])

        #self.plot()
        print("\n")
        print("===================")
        print("END READING 2D")
        print("===================")
        print("\n")
        
    def store_2d(self):
        """
        Method to store the profiles in UFILE format
        """
        try:
            self.rsig
        except:
            self.get_signal()
        for k in self.signals.keys():
            uf_d = {'pre': self.signals[k]['prefix'],
                    'ext': self.signals[k]['suff'],
                    'shot': self.shot,
                    'grid': {'X': {'lbl': tlbl,
                                   'arr': self.rsig[k]['time']},
                             'Y': {'lbl': ylbl,
                                   'arr': self.rsig[k]['rho']}},
                    'data': {'arr': self.rsig[k]['signal'],
                             'lbl': self.signals[k]['dlbl']}}
            ufiles.WU(uf_d, udir=self.path)

    def plot(self):
        for i in self.signals.keys():
            data = self.rsig[i]['signal']
            self._plot_input(data, i)
        

    def _plot_input(self, *argv):
        """
        Function that plots 2D profiles. 
        Useful for benchmark of input quantities

        """
#        if __name__ == '__main__':
#            self.pltframe = tk.TK()
#        else:
#            self.pltframe = tk.Toplevel()
#
#        self.pltframe.title(str(argv[-1]))

        
        n_sp=self.time_u.size
        if n_sp > 10:
            plotind=np.linspace(0, n_sp-1, num=10, dtype=int)
        else:
            plotind = range(n_sp)

        fig = plt.figure(figsize=(5,5), dpi=100)
        ax = fig.add_subplot(111)
        ax.set_title(argv[1])
        for i in plotind:
            #print(i, n_sp)
            #try:
                ax.plot(self.rho, argv[0][i,:], label=self.time_u[i]) 
            #except:
            #    ax.plot(np.linspace(0, 1, num=self.n_rho_eq), \
            #            argv[0][i,:],label=self.time_u[i]) 
        #fig.tight_layout()
        plt.legend(loc='best')
        plt.show()
        #fig=plt.gcf()
#        canvas = FigureCanvasTkAgg(fig, self.pltframe)
#        canvas.show()
#        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#        button = sty.myButton(self.pltframe, 'Dismiss', self.pltframe.destroy, bgc=tksty.qtcol)
#        button.pack()
