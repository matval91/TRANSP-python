import utils.plot_utils as au
import matplotlib.pyplot as plt
import pytransp.trutils.transp_utils as tu
import pytransp.classes.transp_exp as te
import numpy as np

col, col2, styles, my_cmap, dpi = au.define_colors()

def summary_plot(fname_sim=['/home/vallar/TCV/58823/58823V71.CDF'], t=0.):
    """
    Summary plot of a TRANSP simulation
    
    Params:
        fname   str : CDF filename of the simulation

    """
    if type(fname_sim)==str: fname_sim=[fname_sim]

    ncol=4; nrow=3
    f, ax = plt.subplots(nrow, ncol, figsize=[4*ncol,3*nrow])

    for i, fname in enumerate(fname_sim):
        try:
            if fname[-4:]!='.CDF': fname+='.CDF'
            o=te.transp_exp(fname)
        except:
            print('No fname '+fname)
            continue
        if t==0.:
            ind = int(np.size(o.t)/2)
        else:
            ind = np.argmin(o.t-t<0.)

        ax[0,0].plot(o.t, o.ip*1e-3,  col[i], lw=2.3, label=fname[-12:-4]); 
        ax[0,0].set_xlabel(r't (s)'); ax[0,0].set_ylabel(r'$I_{TOT} (kA)$')
        ax[0,0].legend(loc='best')
        ax[0,0].axvline(x=o.t[ind])

        ax[1,0].plot(o.rho[ind,:], o.file['Q'][ind,:], col[i], lw=2.3, label='t='+str(round(o.t[ind],3))+' s');
        if i==0:
            ax[1,0].plot([0., 1.], [1.,1.], 'k--', lw=2.3)
        ax[1,0].set_xlabel(r'$\rho$'); ax[1,0].set_ylabel(r'q')
        ax[1,0].legend(loc='best')
        
        ax[0,1].plot(o.rho[ind,:], o.kin_vars['ne'][ind,:], col[i], lw=2.3);
        ax[0,1].set_xlabel(r'$\rho$'); ax[0,1].set_ylabel(r'$n_e (m^{-3})$')    

        ax[1,1].plot(o.rho[ind,:], o.kin_vars['te'][ind,:]*1e-3, col[i], lw=2.3);
        ax[1,1].set_xlabel(r'$\rho$'); ax[1,1].set_ylabel(r'$T_e (keV)$') 

        ax[0,2].plot(o.rho[ind,:], o.nb_FKP_vars['n'][ind,:], col[i], lw=2.3);
        ax[0,2].set_xlabel(r'$\rho$'); ax[0,2].set_ylabel(r'$n_{fast} (m^{-3})$') 

        ax[1,2].plot(o.rho[ind,:], o.nb_FKP_vars['nbcd'][ind,:]*1e-3, col[i], lw=2.3);
        ax[1,2].set_xlabel(r'$\rho$'); ax[1,2].set_ylabel(r'$j_{fast} (kA\cdot m^{-2})$') 
    
        ax[2,0].plot(o.t, o.file['YAXIS'][:], col[i], lw=2.3)
        ax[2,0].set_xlabel(r't (s)'); ax[2,0].set_ylabel(r'$z_{axis} (cm)$') 
        ax[2,0].axvline(x=o.t[ind])

        ax[2,1].plot(o.rho[ind,:], o.file['TAUPI'][ind,:], col[i], lw=2.3)
        ax[2,1].set_xlabel(r'$\rho$'); ax[2,1].set_ylabel(r'$\tau^{conf}_{ion} (s)$') 
        
        o._calculate_n0()
        ax[2,2].semilogy(o.rho[ind,:], o.n0_tot[ind,:], col[i], lw=2.3)
        ax[2,2].set_xlabel(r'$\rho$'); ax[2,2].set_ylabel(r'$n_0 (m^{-3})$')        
        
        ax[0,3].plot(o.t, o.nb_in_vars['P_D']*1e-6, col[i],marker='x', lw=2.3)
        ax[0,3].set_xlabel(r't (s)'); ax[0,3].set_ylabel(r'$P_{inj} (MW)$')        
        ax[0,3].axvline(x=o.t[ind])

        ax[1,3].plot(o.t, o.nb_in_vars['E_D']*1e-3, col[i], marker='x', lw=2.3)
        ax[1,3].set_xlabel(r't (s)'); ax[1,3].set_ylabel(r'$E_{inj} (keV)$')        
        ax[1,3].axvline(x=o.t[ind])
        
        o._calculate_wmhd()
        ax[2,3].plot(o.t, o.wtot*1e-3, col[i], lw=2.3)
        ax[2,3].plot(o.t, o.wth*1e-3, col[i], ls='--', lw=2.3)        
        ax[2,3].set_xlabel(r't (s)'); ax[2,3].set_ylabel(r'$W (kJ)$')        
        ax[2,3].axvline(x=o.t[ind])        
    for el in ax: 
        for el2 in el:
            el2.grid('on')
    f.tight_layout()
    plt.show()
    return f,ax
