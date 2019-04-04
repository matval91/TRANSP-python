import transp_heating as th
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import define_colors
colors,_,_,_,_ = define_colors()

def compare_nb(fnames=['/home/vallar/TCV/58823/58823V68.CDF','/home/vallar/TCV/58823/58823V69.CDF'],\
                labels=['',''], time=1.25):

    f=plt.figure(); ax=f.add_subplot(111)
    fn=plt.figure(); axn=fn.add_subplot(111)
    fi=plt.figure(); axi=fi.add_subplot(111)

    if labels[0]=='':
        labels = [fname[-12:-4] for fname in fnames]
    for i, fname in enumerate(fnames):
        o=th.transp_heating(fname)
        if o.flag_beam!=1:
            print('No beam in simulation ', labels[i])
            continue
        ind = np.argmin(o.t - time < 0.)
        rho = o.rho[ind,:]; pbe = o.nb_FKP_vars['pe'][ind,:]; n = o.nb_FKP_vars['n'][ind,:]
        pbi = o.nb_FKP_vars['pi'][ind,:];    
        ax.plot(rho,  pbe, colors[i], label=labels[i])
        axi.plot(rho, pbi, colors[i], label=labels[i])
        axn.plot(rho, n*1e6, colors[i], label=labels[i])

    ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'pbe')
    ax.legend(loc='best')
    f.tight_layout(); ax.grid('on');

    axn.set_xlabel(r'$\rho$'); axn.set_ylabel(r'n [1/m$^3$]')
    axn.legend(loc='best')
    fn.tight_layout(); axn.grid('on');

    axi.set_xlabel(r'$\rho$'); axi.set_ylabel(r'pbi')
    axi.legend(loc='best')
    fi.tight_layout(); axi.grid('on');

    plt.show()
