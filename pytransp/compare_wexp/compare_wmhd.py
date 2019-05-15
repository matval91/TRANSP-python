import pytransp.classes.transp_exp as te
import numpy as np
import matplotlib.pyplot as plt
import utils.plot_utils as pu

colors,_,_,_,_ = pu.define_colors()

def compare_wmhd(fname_sim=['/home/vallar/TCV/58823/58823V71.CDF'], fname_exp_w='/home/vallar/TCV/58823/58823wmhd.dat', offset=0.):
    """
    Compares energy from TRANSP with experimental energy

    Parameters:
        fname_sim    arr  : array of strings with the name of CDF files
        fname_exp_w  str  : string of name where energy is stored
        offset       float: offset to exp energy: wmhd=wmhd-offset
    Returns:
        wmhd_exp  dict: dict with time and wmhd of experiment
        wmhd_sim  dict: dict with time and wmhd of CDF files
    """
    if type(fname_sim)==str: fname_sim=[fname_sim]

    wmhddata = np.loadtxt(fname_exp_w)
    
    f=plt.figure()
    ax=f.add_subplot(111)
    ax.plot(wmhddata[:,0], wmhddata[:,1]*1e-3-offset, 'm--', label='EXP')
    wmhd_exp = {}; wmhd_exp['t']=wmhddata[:,0]; wmhd_exp['wmhd']=wmhddata[:,1]
    wmhd_sim = {};
    for i, fname in enumerate(fname_sim):
        try:
            if fname[-4:]!='.CDF': fname+='.CDF'
            o=te.transp_exp(fname)
        except:
            print('No file named ', fname)
            continue
        o._calculate_wmhd()
        ax.plot(o.t, o.wtot*1e-3, colors[i], label=fname[-12:-4])
        ax.plot(o.t, o.wth*1e-3, colors[i], linestyle='--', label=fname[-12:-4])
        wmhd_sim[fname[-12:-4]] = {}
        wmhd_sim[fname[-12:-4]]['t']    = o.t
        wmhd_sim[fname[-12:-4]]['wmhd'] = o.wtot
        wmhd_sim[fname[-12:-4]]['wth']  = o.wth
    ax.set_xlabel(r'Time (s)'); ax.set_ylabel(r'W$_{MHD}$ [kJ]')
    ax.legend(loc='best')
    f.tight_layout()
    ax.grid('on')
    plt.show()
    return wmhd_exp, wmhd_sim

