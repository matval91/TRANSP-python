import utils.plot_utils as au
import matplotlib
import matplotlib.pyplot as plt
import pytransp.trutils.transp_utils as tu
col, col2, styles, my_cmap, dpi = au.define_colors()

def plot_eq(to, time=[0], f=0):
    """ plot equilibrium

    Plots equilibrium

    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        None    
    Note:
    """        
    au.common_style()
    ind = tu._time_to_ind(to.t, time)
    ind=ind[0]

    flag_leg=1
    ls='-'
    if not isinstance(f, matplotlib.figure.Figure):
        f = plt.figure(figsize=(12,10))
        axj = f.add_subplot(221)
        axp = f.add_subplot(222)
        axf = f.add_subplot(223)
    else:
        axj, axp, axf, axq=f.axes
        flag_leg=0
        ls='--'
    axj.plot(to.rho[ind,:], to.eq_vars['j'][ind,:]*1e-3, 'k', label=r'$j_{TOT}$', lw=2., linestyle=ls)
    axj.plot(to.rho[ind,:], to.eq_vars['joh'][ind,:]*1e-3, 'b', label=r'$j_{OH}$', lw=2., linestyle=ls)
    axj.plot(to.rho[ind,:], to.eq_vars['jbs'][ind,:]*1e-3, 'r', label=r'$j_{BS}$', lw=2., linestyle=ls)
    add_curr=to.eq_vars['j'][ind,:]-to.eq_vars['joh'][ind,:]-to.eq_vars['jbs'][ind,:]
    axj.plot(to.rho[ind,:], add_curr*1e-3, 'g', label=r'$j_{Add. heat.}$', lw=2., linestyle=ls)

    axp.plot(to.rho[ind,:], to.eq_vars['p'][ind,:]*1e-3, 'k', label=r'Total', lw=2., linestyle=ls)
    axp.plot(to.rho[ind,:], to.eq_vars['pth'][ind,:]*1e-3, 'r', label=r'Th.', lw=2., linestyle=ls)
    axp.plot(to.rho[ind,:], to.eq_vars['pnth'][ind,:]*1e-3, 'b', label=r'Non th.', lw=2., linestyle=ls)
    
    axf.plot(to.rho[ind,:], to.eq_vars['pol_flux'][ind,:]*10., 'k', label=r'Pol. flux x 10', lw=2., linestyle=ls)
    axf.plot(to.rho[ind,:], to.eq_vars['tor_flux'][ind,:], 'b', label=r'Tor. Flux', lw=2., linestyle=ls)
    axq = axf.twinx()  # instantiate a second axes that shares the same x-axis
    axq.plot(to.rho[ind,:], to.eq_vars['q'][ind,:], 'r', label=r'q', lw=2., linestyle=ls)
    axq.plot([0,1], [1,1], 'r--')
    axq.set_ylim([0,10]); axf.set_ylim([0,0.5])

    #========================================================
    # SET TICK LOCATION
    #========================================================
    if flag_leg==1:
        axq.tick_params(axis='y', labelcolor='r')
        axq.set_ylabel(r'q', color='r')
        au.limit_labels(axj, r'$\rho_{TOR}$', r'j [kA/m$^2$]','' )
        au.limit_labels(axp, r'$\rho_{TOR}$', r'p [kPa]','' )
        axf.set_xlabel(r'$\rho_{TOR}$')
        axf.set_ylabel(r'Fluxes (Wb/rad)')
        axj.legend(loc='best')
        axf.legend(loc='upper left'); axf.grid('on')
        axp.legend(loc='upper right'); axp.grid('on')
        f.tight_layout()
        f.subplots_adjust(top=0.9)
    else:
        axq.set_yticklabels([])
    try:
        titles=f.get_children()[-1]
        title=titles.get_text()
        if titles.get_text() == '':
            f.suptitle('{} t={:.2f} s'.format(to.fname[-12:-4], time))
        else:
            newtitle = title+' {}'.format(to.fname[-12:-4])
            titles.set_text(newtitle)
    except:
        f.suptitle('{} t={:.2f} s'.format(to.fname[-12:-4], time))

    plt.show()    
