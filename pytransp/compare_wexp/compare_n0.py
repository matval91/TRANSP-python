import pytransp.classes.transp_exp as te
import numpy as np
import matplotlib.pyplot as plt
import utils.plot_utils as pu

colors,_,_,_,_ = pu.define_colors()

def compare_n0(fname_sim=['/home/vallar/TCV/58823/58823V71.CDF'], fname_exp='/home/vallar/TCV/58823/58823n0.dat', t=0.9):
    """
    Compares edge density between CDF output of TRANSP and n0 from experiments
    
    Parameters:
        fname_sim  arr  : array of strings with the name of CDF files
        fname_exp  str  : string of name where n0 at edge is stored
        t          float: time where to compute profile
    Returns:
        n0_exp  dict: dict with time and n0 of experiment
        n0_sim  dict: dict with time and n0 of CDF files
    
    """
    if type(fname_sim)==str: fname_sim=[fname_sim]

    # read data from experiments
    n0data = np.loadtxt(fname_exp)

    # start plot, it will be needed inside the for loop
    f1d = plt.figure(); ax1d=f1d.add_subplot(111)
    f=plt.figure(); ax=f.add_subplot(111) 
    ax.plot(n0data[:,0], n0data[:,1], 'k--', label='EXP')
    n0_exp = {}; n0_exp['t']=n0data[:,0]; n0_exp['n0']=n0data[:,1]
    n0_sim = {};   
    #plots for every filename chosen
    for i, fname in enumerate(fname_sim):
        try:
            if fname[-4:]!='.CDF': fname+='.CDF'
            o=te.transp_exp(fname)
        except:
            continue
        ind = np.argmin(o.t-t < 0.)
        o._calculate_n0()
        n0_sim[fname[-12:-4]] = {}
        n0_sim[fname[-12:-4]]['t']  = o.t
        n0_sim[fname[-12:-4]]['n0'] = o.n0_tot[:,-1]
        ax.plot(o.t, o.n0_tot[:,-1]*1e6, colors[i], label=fname[-12:-4])
        ax1d.semilogy(o.rho[ind,:], o.n0_tot[ind,:]*1e6, colors[i], label=fname[-12:-4])
        #ax1d.semilogy([1., 1.1], [2e16, 2e16], colors[i])

    #making nice plots
    ax.set_xlabel(r'Time (s)'); ax.set_ylabel(r'n$_0^{edge}$ [1/m$^3$]')
    ax.legend(loc='best')
    ax.set_ylim([0, 2e16]); ax.set_xlim([0, 2.])
    f.tight_layout(); ax.grid('on');

    ax1d.set_xlabel(r'$\rho$'); ax1d.set_ylabel(r'n$_0$ [1/m$^3$]')
    ax1d.legend(loc='best')
    f1d.tight_layout(); ax1d.grid('on')
    plt.show()
    return n0_exp, n0_sim


