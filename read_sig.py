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
        self._univec_channel = {'EC':{},\
                                'NBH':{}, \
                                'DNB':{}}  # this is a dictionary of object
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
#                                'string': r'tcv_bphi()',
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
        self.chann_sig = {'EC': {'string': r'\results::toray_input.pgyro',
                                 'yl': r'P$_{EC}$  [W]',
                                 'ytk': [0, 1e6, 2e6, 3e6],
                                 'chlab': ['L1','L2','L3','L4','L5',\
                                        'L6','L7','L8','L9'],
                                 'lbl': 'EC power'.ljust(20) + 'W',
                                 'suff': 'EC', 'prefix': 'P'},
                          'NBH': {'string': r'\atlas::nbh.data.main_adc:data',
                                  'yl': r'P$_{NBI}$  [W]',
                                  'ytk': [0, 4e5, 8e5, 1e6],
                                  'chlab': ['NBH'],
                                  'lbl': 'NBI power'.ljust(20) + 'W',
                                  'suff': 'NBI', 'prefix': 'P'},
                          'DNB': {'string': r'tofill',
                                  'yl': r'P$_{DNB}$  [W]',
                                  'ytk': [0, 4e5, 8e5, 1e6],
                                  'chlab': ['DNB'],
                                  'lbl': 'DNB power'.ljust(20) + 'W',
                                  'suff': 'NBI', 'prefix': 'P'}}


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
#            if stri == 'tcv_bphi()':
#                rax = self.conn.tdi(r'\results::r_axis')
#                rax = rax.values[_iidx]
#                dat *= rax
            tim = tim[_iidx]
            dummy = interpolate.InterpolatedUnivariateSpline(tim, \
                                                             dat, ext=0)
            self._univec[k] = dict([('spline', dummy)])
            
        #self._read_channels()
        # add also the univariate spline interpolation for the NBI
        print("\n")
        print("===================")
        print("END READING 1D")
        print("===================")
        print("\n")

        
        
    def _read_channels(self):
        """
        Reads quantities 1D but with channels (EC, [NBI+DNBI])
        """

        #EC
        print('Reading signal '+self.chann_sig['EC']['string'])
        x,y = self._read_ecnbi('EC')
        for i, channels in enumerate(self.indGyro):
            ll=self.chann_sig['EC']['chlab'][channels]
            dummy=interpolate.InterpolatedUnivariateSpline(x, y[:,i], ext=0)
            self._univec_channel['EC'][ll] = dict([('spline', dummy)])
           
        #NBH... to add DNBI
        print('Reading signal '+self.chann_sig['NBH']['string'])
        x,y = self._read_ecnbi('NBH')
        dummy = interpolate.InterpolatedUnivariateSpline(x, y, ext=0)
        self._univec_channel['NBH'] = dict([('spline', dummy)])

        #DNBI
        print('Reading signal '+self.chann_sig['DNB']['string'])
        x,y = self._read_ecnbi('DNB')
        dummy = interpolate.InterpolatedUnivariateSpline(x, y, ext=0)
        self._univec_channel['DNB'] = dict([('spline', dummy)])        
        
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
        #self.rsig_channel = {'EC':{},'NBH':{}, 'DNB':{}}
        self.rsig_channel = {'NBH':{}, 'DNB':{}}

        for k in self._univec.keys():
            self.rsig[k] = dict([('data',
                                  self._univec[k]['spline'](self.time_u)),
                                 ('time', self.time_u)])
                            
        # for k in self.rsig_channel.keys():
        #     self.rsig_channel[k] = dict([
        #             ('data', self._univec_channel[k]['spline'](self.time_u)),\
        #             ('time', self.time_u)
        #             ])  
            


    def set_zeff(self):
        """
        Allows to set a constant value of Zeff
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
        self._store_1d_channels()
            
        
    def _store_1d(self):
        """
        This method store the data which has already been read. It does not
        store the ECRH which have multiple sources
        """

        # first of all check the existence of self. rsig
        try:
            self.rsig
        except:
            self._read_1d()

        # now we cycle to save the appropriate ufiles
        # MV 10/2016: for P<shot>.NBI you must write also the number of NB lines
        # used, so it has a separate writing. This is due to nubeam_driver
        for k in self.signals.keys():
            uf_d = {'pre': self.signals[k]['prefix'],
                    'ext': self.signals[k]['suff'],
                    'shot': self.shot,
                    'grid': {'X': {'lbl': tlbl,
                                   'arr': self.rsig[k]['time']}},
                    'data': {'lbl': self.signals[k]['lbl'],
                             'arr': self.rsig[k]['data']}}
            
            ufiles.WU(uf_d, udir=self.path)
        
        
    def _store_1d_channels(self):
        """
        Stores NBI and EC
        """
        self._store_nbi()
        if 'EC' in self.chann_sig.keys():
            self._store_ech()
        
    def read_1d(self):
        """
        Returns the dictionary of the signals (Ip,RBT,ZEFF,VLOPP,NBI)
        """
        try:
            self.rsig
        except:
            self._read_1d()
        #self.plot_input()
        #return self.rsig

    def _read_ecnbi(self, label):
        """
        Makes possible to read or EC or NBI, which are 1D signals but must be 
        written to uFiles as 2D signals (time, channel)
        """
        if label == 'NBH':
            x,y=self._read_nbi()
        elif label=='DNB':
            x,y = self._read_dnb()
        else:
            x,y=self._read_ecrh()
        
        return x,y
            
    def _get_nbi(self):

        """
        Get and save the corresponding ufile for the NBI signal.
        Mimic the Matlab script from Alexander Karpushov

        Parameters
        ----------
        shot

        Returns
        ----------
        Returns the xdata for the NBI as read in \atlas

        """

        # read the estimated neutral power
        # the data on ATLAS are written as MW, thus the conversion factor
        # (*1e6) has been added (M. Vallar 10/2016)
        self.nbiC = self.conn.tdi(r'\atlas::nbh.data.main_adc:data')*1e6

    def _read_nbi(self):

        """
        For the NBI we get as an output directly the time and the power

        """
        try:
            self.nbiC
        except:
            self._get_nbi()

        # limit our self in the time interval chosen
        _iidx = ((self.nbiC.dim_0.values >= self.tbeg-self.dt) &
                 (self.nbiC.dim_0.values <= self.tend+self.dt))

        # this is the power of the neutrals, i.e. ion power x neutral. efficiency
        return self.nbiC.dim_0.values[_iidx], self.nbiC.values[_iidx, 36]
 
    def _read_dnb(self):

        """
        For the DNB we get as an output directly the time and the power
        TO FIX BETTER; NOW INITALISES TO 0
        """
#        try:
#            self.nbiC
#        except:
#            self._get_dnb()

#        # limit our self in the time interval chosen
#        _iidx = ((self.nbiC.dim_0.values >= self.tbeg-self.dt) &
#                 (self.nbiC.dim_0.values <= self.tend+self.dt))
#
#        # this is the power of the neutrals, i.e. ion power x neutral. efficiency
#        return self.nbiC.dim_0.values[_iidx], self.nbiC.values[_iidx, 36]
        time_dnb  = np.linspace(self.tbeg, self.tend, num=100, dtype=float)
        power_dnb = np.zeros(100, dtype=float)
        return time_dnb, power_dnb
                              
    def _neuteff(self,E):
        """
        Calculates neutralization efficiency as 
        polynomial fit of degree 5 of ELPRO data
                               AK (2003)
        """
        
        P = [-1.3785841564824268e-10, \
             2.48797392816381106e-08,\
             -1.151346804847086e-07,\
             -0.00015463763952549783,\
             -0.0015147724400023021,\
             0.89390606995761313]

        eff = np.polyval(P,E);
        return eff
    
    def _nbi_energy(self):
        """
        Function to get average energy and energy fraction of the NBI
        """
        try:
            self.nbiC
        except:
            self._get_nbi()
        #power_mix = [0.773  ,0.162, 0.063, 0.002]

        einj_arr=self.nbiC[:,32].data
        ptot = self.nbiC[:,36]
        index = ptot > 0.5e6
        self.einj_NBI = np.mean(einj_arr[index])/1.e6 #energy in keV
        self.effNeu = self._neuteff(self.einj_NBI)
    
    def _store_nbi(self):
        """
        Store the appropriate U-files for the NBI
        """
        try:
            self.nbiC
        except:
            self._get_nbi()
            self._read_1d()
            
            
        NBIdict=self.chann_sig['NBH']
        data=np.array([self.rsig_channel['NBH']['data'], \
              self.rsig_channel['DNB']['data']]).T
        time=self.rsig_channel['NBH']['time']

        arr_nbi=np.array([1.0, 2.0])
        xlbl = 'Channel number'
        dlbl = NBIdict['lbl']
        pref = NBIdict['prefix']
        suff = NBIdict['suff']
        uf_d = {'pre': pref, 'ext': suff, 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl,
                               'arr': time},
                         'Y': {'lbl': xlbl,
                               'arr': arr_nbi}},
                'data': {'lbl': dlbl, 'arr': data}}
        ufiles.WU(uf_d, udir=self.path)

    def _store_ech(self):
        """
        Store the appropriate U-files for the ECH
        """
        try:
            self.echC
        except:
            self.read_ecrh()
            
        ECdict = self.chann_sig['EC']
        try:
            time = self.rsig_channel['EC']['L1']['time']
        except:
            time = self.rsig_channel['EC']['L4']['time']
            
        data = np.zeros((len(time), len(self.indGyro)), dtype=float)
        for i,el in enumerate(self.rsig_channel['EC'].keys()):
            data[:, i] = self.rsig_channel['EC'][el]['data']
            data[data[:,i]<1.,i]=0.

        arr_ec = np.arange(len(self.indGyro), dtype=int)+1
        xlbl = 'Channel number'
        dlbl = ECdict['lbl']
        pref = ECdict['prefix']
        suff = ECdict['suff']
        uf_d = {'pre': pref, 'ext': suff, 'shot': self.shot,
                'grid': {'X': {'lbl': tlbl,
                               'arr': time},
                         'Y': {'lbl': xlbl,
                               'arr': arr_ec}},
                'data': {'lbl': dlbl, 'arr': data}}
        ufiles.WU(uf_d, udir=self.path)
        
    def _get_ecrh(self):

        """
        Hidden method for reading the gyrotron power delivered

        """
        #self.echC = self.conn.tdi(r'\results::toray.input:p_gyro')
        self.echC = self.tree.getNode(r'\results::toray.input:p_gyro')

    def _read_ecrh(self):
        """
        Return the power of the gyrotrons in the choosen time interval
        """
        try:
            self.echC
        except:
            self._get_ecrh()

        #_iidx = ((self.echC.dim_0.values >= self.tbeg-self.dt) &
        #         (self.echC.dim_0.values <= self.tend+self.dt))
#        _iidx = np.argwhere(((self.echC.getDimensionAt(0).value >= self.tbeg-self.dt) &
#                 (self.echC.getDimensionAt(0).value <= self.tend+self.dt)))
        _indGyro = np.argwhere(~np.isnan(np.nanmean(self.echC.data(), axis=1)))
        self.indGyro = _indGyro[:-1,0]
#        
#        return self.echC.getDimensionAt(0).value[_iidx], \
#            np.nan_to_num(self.echC.data()[self.indGyro,_iidx])
        time_dnb  = np.linspace(self.tbeg, self.tend, num=100, dtype=float)
        power_dnb = np.zeros(100, dtype=float)
        return time_dnb, power_dnb
    
    
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
        self.time_u = np.arange(self.tbeg, self.tend, self.dt)
        if self.tbeg==self.tend:
            self.time_u = np.array([self.tbeg])
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
        self.signals = {#'ne': {'string': r'\results::conf:ne',
                        'ne': {'string': r'\tcv_shot::top.results.thomson.profiles.auto:ne',
                               'zl': r'n$_e$ [10$^{19}$m$^{-3}$]',
                               'yl': r'$\rho_{\phi}$',
                               'suff': 'ELE', 'prefix': 'N',
                               'xlbl': tlbl, 'ylbl': ylbl,
                               'dlbl': 'NE                  [cm^-3]  '},
                        #'te': {'string': r'\results::conf:te',
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
                        #'vtor': {'string': r'\results::cxrs.proffit:vi_tor',
                        #         'zl': r'v$_{\phi}$ [rad/s]',
                        #         'yl': r'$\rho_{\phi}$',
                        #         'suff': 'TOR', 'prefix': 'V',
                        #         'xlbl': tlbl,
                        #         'ylbl': ylbl,
                        #         'dlbl': 'omg'.ljust(23) + '[rad/s]'},
                        #'p': {'string': r'\results::conf:pe',
#                        'p': {'string':r'\tcv_shot::top.results.thomson.profiles.auto:pe',
#                                 'zl': r'P [Pa]',
#                                 'yl': r'$\rho_{\phi}$',
#                                 'suff': 'EQ', 'prefix': 'P',
#                                 'xlbl': tlbl,
#                                 'ylbl': ylbl,
#                                 'dlbl': 'P'.ljust(23) + '[Pa]'}}
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
            # retrieve the time basis of the data CONF
            if k == 'ti':
                tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(1).data()
            else:
                #this will work with THOMSON:
                tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()
           
            # choose only the timining between tbeg and tend
            if k=='vtor':
                _idx=((tim >= self.tbeg-10*self.dt) & (tim <= self.tend+10*self.dt))
            else:
                _idx = ((tim >= self.tbeg) & (tim <= self.tend))
            if len(tim[_idx])<=4:
                _idx = ((tim >= self.tbeg-0.1) & (tim <= self.tend+0.1))
            #print("TIM LEN", len(tim[_idx]))
            #print(_idx)
            #data = self.tree.getNode(self.signals[k]['string']).data()[_idx, :]
            if k!='ti':
            # with thomson
                data = self.tree.getNode(self.signals[k]['string']).data()[:, _idx]
            else:
                data = self.tree.getNode(self.signals[k]['string']).data()[_idx,:]
            #data = np.mean(data, axis=0)
            if k=='ne':
                data=data*1e-6
            # and now limit the time
            tim = tim[_idx]
            #print("TIM", tim)
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
            #if k == 'per':
                rhop = self.tree.getNode(r'\results::thomson.profiles.auto:rho').data()
                rhot = np.zeros((tim.size, rhop.shape[0]))
                for t, i in zip(tim, range(tim.size)):
                    rhot[i, :] = self.eq.psinorm2phinorm(rhop, t,
                                                         sqrt=True)
            #else:
            #    # load the rho toroidal saved in the conf tree
            #    rhot = self.tree.getNode(r'\results::conf:rhotor').data()[_idx, :]

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
