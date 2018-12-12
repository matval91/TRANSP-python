import a4py.utils.ascot_utils as au
import matplotlib.pyplot as plt
import pytransp.utils.transp_utils as tu
col, col2, styles, my_cmap, dpi = au.define_colors()

def plot_input(to, time):
    """ plot all the input

    Plots all the input parameters (vol-avgd and profiles)
    calling the two relative functions

    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        None    
    Note:
    """
    #ind = self._time_to_ind(time)
    plot_input_1d(to, time)
    plot_input_prof(to, time)


def plot_input_1d(to, time, axne, axTe, ls):
    """ plot 1d input

    Plots ne,ni,te,ti vol-avgd

    Parameters:
        time (arr) : array with the times where to plot the lines
    Attributes:
        None    
    Note:
    """
    au.common_style()

    if axne==0:
        f = plt.figure(figsize=(8, 8))
        axne = f.add_subplot(211)
        axTe = f.add_subplot(212, sharex=axne)
        fig_flag = 0 # this means it wasn't initialized 
    #axzf = f.add_subplot(313, sharex=axne)
    else:
        f=plt.gcf()
        fig_flag=1

    t=to.t
    for y, c, l in zip([to.ne_vavg, to.ni_vavg], ['k','r'], [r'e', r'i']):
        au._plot_1d(t, y, ax=axne, color=c, label=l, ls=ls)
    for y, c, l in zip([to.Te_vavg, to.Ti_vavg], ['k','r'], [r'e', r'i']):
        au._plot_1d(t, y, ax=axTe, color=c, label=l, ls=ls)

    if fig_flag==0:
        #====================================================
        # Correct ticks and xylabels
        #====================================================
        au.limit_labels(axne, r'Time [s]', r'$\langle n \rangle$ [$1/m^3$]','' )
        au.limit_labels(axTe, r'Time [s]', r'$\langle T \rangle$ [$eV$]','' )

    #====================================================
    # Plot vertical lines
    #====================================================
    ind = tu._time_to_ind(to.t, time)
    for ax in [axne, axTe]:#, axzf]:
        if len(ind)!=1:
            for i, el in enumerate(ind):
                ax.axvline(x=to.t[el], color=col[i], lw=2., linestyle='--')
        if fig_flag==0:
            ax.legend(loc='best')

    f.tight_layout()
    plt.show()    


def plot_input_prof(to, time=[0]):
    """ plot 1d input

    Plots input quantities, such as ne, ni, Te, Ti profiles

    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        None    
    Note:
    """        
    au.common_style()

    f = plt.figure(figsize=(12,10))
    axne = f.add_subplot(221)
    axni = f.add_subplot(222, sharey=axne)
    axTe = f.add_subplot(223)
    axTi = f.add_subplot(224, sharey=axTe)
    ind = tu._time_to_ind(to, time)

    for i, el in enumerate(ind):
        lab = r't = {:.2f}'.format(to.t[el])
        axne.plot(to.rho[el,:], to.kin_vars['ne'][el,:], col[i], label=lab, lw=2.)
        axni.plot(to.rho[el,:], to.kin_vars['ni'][el,:], col[i], label=lab, lw=2.)
        axTe.plot(to.rho[el,:], to.kin_vars['te'][el,:]*1e-3, col[i], label=lab, lw=2.)
        axTi.plot(to.rho[el,:], to.kin_vars['ti'][el,:]*1e-3, col[i], label=lab, lw =2.)

    #========================================================
    # SET TICK LOCATION
    #========================================================
    au.limit_labels(axne, r'$\rho$', r'$n_e$ [$1/m^3$]','' )
    au.limit_labels(axTe, r'$\rho$', r'$T_e$ [$keV$]','' )
    au.limit_labels(axni, r'$\rho$', r'$n_i$ [$1/m^3$]','' )
    au.limit_labels(axTi, r'$\rho$', r'$T_i$ [$keV$]','' )

    axTi.legend(loc='best')

    f.tight_layout()
    #f.subplots_adjust(left=0.2)
    plt.show()    
