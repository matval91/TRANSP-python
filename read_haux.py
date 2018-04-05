from __future__ import print_function
import os
import ufiles
import tcv
import MDSplus as mds
import numpy as np
from scipy import interpolate
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
class AUX:
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

        self.time_u = np.arange(self.tbeg, self.tend, self.dt)

        self.chann_sig = {'EC': {'string': r"\tcv_shot::top.results.toray.input:p_gyro",
                                 'yl': r'P$_{EC}$  [W]',
                                 'ytk': [0, 1e6, 2e6, 3e6],
                                 'chlab': ['L1','L2','L3','L4','L5',\
                                        'L6','L7','L8','L9'],
                                 'lbl': 'EC power'.ljust(20) + 'W',
                                 'suff': 'EC', 'prefix': 'P'},
                          'NBH': {'string':r"\RESULTS::NBH:POWR_TCV",
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
        stri_launch = r'\tcv_shot::top.ecrh.measurements.launchers'
        self.X3_angles = {
            'L7':{'theta':
                  {'string':stri_launch+'.theta:x3_7'},
                  },
            'L8':{'theta':
                  {'string':stri_launch+'.theta:x3_8'},
                  },
            'L9':{'theta':
                  {'string':stri_launch+'.theta:x3_9'},
                  }
            }
                  
        self.X2_angles = {
            'L1':{'theta':
                  {'string':stri_launch+'.theta:x2_1'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_1'},
                  },
            'L2':{'theta':
                  {'string':stri_launch+'.theta:x2_2'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_2'},
                  },
            'L3':{'theta':
                  {'string':stri_launch+'.theta:x2_3'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_3'},
                  },
            'L4':{'theta':
                  {'string':stri_launch+'.theta:x2_4'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_4'},
                  },
            'L5':{'theta':
                  {'string':stri_launch+'.theta:x2_5'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_5'},
                  },
            'L6':{'theta':
                  {'string':stri_launch+'.theta:x2_6'},
                  'phi':
                  {'string':stri_launch+'.phi:x2_6'},
                  },
            }                      


        print("\n")
        print("===================")
        print("Initialisation of AUX HEATING signals Done")
        print("===================")
        print("\n")

    def read_sig(self):
        try:
            self.rsig['NBH']['data'].mean()
        except:
            self._read_sig()
        
    def _read_sig(self):

        """
        We define the appropriate attribute rsig, which is a dictionary
        with the signals and the time at the resolution needed
        """
        # check we define the appropriate
        # dictionary of UnivariateSpline
        # representation
        try:
            self._univec_channel
        except:
            self._get_spline()

        self.rsig = {'EC':{}, 'NBH':{}, 'DNB':{}}

        for i, ch in enumerate(self.indGyro):
            ll=self.chann_sig['EC']['chlab'][ch]
            self.rsig['EC'][ll] = dict([\
                    ('data', self._univec_channel['EC'][ll]['spline'](self.time_u)),\
                    ('time', self.time_u)])
            self.rsig['EC'][ll]['data'] = np.nan_to_num(self.rsig['EC'][ll]['data'])
        for nn in ['NBH','DNB']:
                self.rsig[nn] = dict([\
                     ('data', self._univec_channel[nn]['spline'](self.time_u)),\
                     ('time', self.time_u)])
                
    def _get_spline(self):
        """
        Reads quantities 1D but with channels (EC, [NBI+DNBI])
        """

        self._univec_channel = {'EC':{},\
                                'NBH':{}, \
                                'DNB':{}}  # this is a dictionary of object
        #EC
        print('Reading signal '+self.chann_sig['EC']['string'])
        x,y = self._read_ecrh()
        for i, channels in enumerate(self.indGyro):
            ll=self.chann_sig['EC']['chlab'][channels]
            dummy=interpolate.InterpolatedUnivariateSpline(x, y[i,:], ext=0)
            self._univec_channel['EC'][ll] = dict([('spline', dummy)])
           
        #NBH... to add DNBI
        print('Reading signal '+self.chann_sig['NBH']['string'])
        x,y = self._read_nbi()
        dummy = interpolate.InterpolatedUnivariateSpline(x, y, ext=0)
        self._univec_channel['NBH'] = dict([('spline', dummy)])

        #DNBI
        print('Reading signal '+self.chann_sig['DNB']['string'])
        x,y = self._read_dnb()
        dummy = interpolate.InterpolatedUnivariateSpline(x, y, ext=0)
        self._univec_channel['DNB'] = dict([('spline', dummy)])        
  
                
    def store_aux(self):
        """
        Stores NBI and EC
        """
        self._store_nbi()
        self._store_ech()
        
    def _read_ecrh(self):
        """
        Return the power of the gyrotrons in the choosen time interval
        """
        stri = self.chann_sig['EC']['string']
        echC = self.tree.getNode(stri)
        tech = echC.getDimensionAt(0).data()
        #Filtering on gyrotrons turned on in this shot
        _indGyro = np.argwhere(~np.isnan(np.nanmean(echC.data(), axis=1)))
        self.indGyro = _indGyro[:-1,0]
        return tech, echC.data()[self.indGyro,:]

    def _read_nbi(self):

        """
        For the NBI we get as an output directly the time and the power

        """
        nbiC = self.tree.getNode(self.chann_sig['NBH']['string'])
        tnbi = nbiC.getDimensionAt(0).data()
        pnbi = nbiC.data()*1e6 #power in MW
        return tnbi, pnbi

        # self.nbiC = self.conn.tdi(r'\atlas::nbh.data.main_adc:data')*1e6
        # # limit our self in the time interval chosen
        # _iidx = ((self.nbiC.dim_0.values >= self.tbeg-self.dt) &
        #          (self.nbiC.dim_0.values <= self.tend+self.dt))

        # # this is the power of the neutrals, i.e. ion power x neutral. efficiency
        # return self.nbiC.dim_0.values[_iidx], self.nbiC.values[_iidx, 36]
 
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
            self.rsig['NBH']['data'].mean()
        except:
            self.read_sig()
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
            self.rsig['NBH']['data'].mean()
        except:
            self._read_sig()
            
            
        NBIdict=self.chann_sig['NBH']
        data=np.array([self.rsig['NBH']['data'], \
              self.rsig['DNB']['data']]).T
        time=self.rsig['NBH']['time']

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
            self.rsig['EC'].mean()
        except:
            self.read_sig()
            
        ECdict = self.chann_sig['EC']
        try:
            time = self.rsig['EC']['L1']['time']
        except:
            time = self.rsig['EC']['L4']['time']
            
        data = np.zeros((len(time), len(self.indGyro)), dtype=float)
        for i,el in enumerate(self.rsig['EC'].keys()):
            data[:, i] = self.rsig['EC'][el]['data']
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
        key_arr = [ 'NBI', 'EC']

        for indi, i in enumerate(key_arr):
            ax = fig.add_subplot(r_nsp, c_nsp,indi+1)
            name = self.chann_sig.keys()[indi]
            for jj in self.rsig[name].keys():
                data = self.rsig[name][jj]['data']
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
