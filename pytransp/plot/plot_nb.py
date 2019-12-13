import numpy as np
import utils.plot_utils as au
import matplotlib.pyplot as plt
import pytransp.trutils.transp_utils as tu

col, col2, styles, my_cmap, dpi = au.define_colors()


def plot_deposition(th, time=[0]):
    """
    """
    ind = tu._time_to_ind(th.t, time, inj_index=th.inj_index)
    plot_deposition_1d(th, ind)
    plot_deposition_prof(th, ind)

def plot_deposition_timeslice(th, timeslice):
    """
    THIS ONE IS UNUSED
    """
    ind=np.argmin(th.t-timeslice<0)
    ind=np.array(ind)
    th.au._plot_deposition_prof(ind=[ind])

def plot_deposition_1d(th, time=[0], axp=0, axf=0, ls='-'):
    """
    """
    try:
        th.gi.mean()
    except:
        th._calculate_gi()
    try:
        th.pe.mean()
    except:
        th._calculate_scalars()

    au.common_style()
    ind = tu._time_to_ind(th.t, time, inj_index=th.inj_index)
    if axp==0:
        f = plt.figure(figsize=(10, 8))
        axp = f.add_subplot(211)
        axf = f.add_subplot(212, sharex=axp)
        fig_flag = 0 # this means it wasn't initialized 
    else:
        f=plt.gcf()
        fig_flag=1            
    t=th.t; ptot=th.pi+th.pe
    for y, c, l in zip([th.pe, th.pi], ['k','r'], [r'e', r'i']):
        au._plot_1d(t, y*1e-3, ax=axp, color=c, label=l, ls=ls)
    for y, c, l in zip([th.pe, th.pi], ['k','r'], [r'e', r'i']):
        au._plot_1d(t, y*100./ptot, ax=axf, color=c, label=l, ls=ls)

    axf.plot(th.t, th.gi_vavg*1e2, 'c', lw=2.3, label=r'G')

    if fig_flag==0:
        #====================================================
        # Correct ticks and xylabels
        #====================================================
        au.limit_labels(axp,r'Time [s]',r'P[kW]','')
        au.limit_labels(axf,r'Time [s]',r'%','')

    #====================================================
    # Plot vertical lines
    #====================================================
    for ax in [axp, axf]:
        if ind[0]!=0:
            for i, el in enumerate(ind):
                ax.axvline(x=th.t[th.inj_index[el]], color=col[i], lw=2., linestyle='--')

    axp.set_ylim([0,160]); axf.set_ylim([0,100.])
    f.tight_layout()
    plt.show()     


def plot_deposition_multishot(ths, time):
    """
    plots deposition for multiple simulations
    """
    au.common_style()
    
    f = plt.figure(figsize=(10,10))
    f.text(0.5, 0.9, r't={:.2f}s'.format(time), fontsize=18)
    axe = f.add_subplot(221)
    axi = f.add_subplot(222, sharex=axe, sharey=axe)
    axn = f.add_subplot(223, sharex=axe)
    axj = f.add_subplot(224, sharex=axe)
    for index,th in enumerate(ths):
        i = tu._time_to_ind(th.t, time, inj_index=th.inj_index)
        i=i[0]
        lab=r'{:s}'.format(th.runid)
        x=th.rho[i,:]
        axi.plot(x, th.nb_FKP_vars['pi'][i,:]*1e-3, col[index], lw=2, label=lab)
        axe.plot(x, th.nb_FKP_vars['pe'][i,:]*1e-3, col[index], lw=2, label=lab)
        axn.plot(x, th.nb_FKP_vars['n'][i, :]/th.kin_vars['ne'][i,:]*100., col[index], lw=2, label=lab)
        axj.plot(x, th.nb_FKP_vars['nbcd'][i,:]*1e-3, col[index],lw=2.,label=lab)        
    au.limit_labels(axe,r'$\rho$', r'$P_e$ [$kW/m^3$]','')
    au.limit_labels(axi,r'$\rho$', r'$P_i$ [$kW/m^3$]','')
    au.limit_labels(axn,r'$\rho$', r'$n_f/n_e$ [%]','')
    au.limit_labels(axj,r'$\rho$', r'j shielded [$kA/m^2$]','')

    axe.legend(loc='best')
    #f.text(0.01, 0.01, th.fname)

    f.tight_layout()
    #f.subplots_adjust(right=0.8)
    plt.show()    



def plot_deposition_prof(th, ind=[0]):
    """
    Plots deposition to ions, electrons
    """
    au.common_style()
    try:
        th.nb_FKP_vars['pi'].mean()
    except:
        print("No data")
        return
    if len(ind)==1 and ind[0]==0:
        ind = np.linspace(0, len(th.inj_index)-1, 5, dtype=int)
    f = plt.figure(figsize=(18,6))
    axe = f.add_subplot(131)
    axi = f.add_subplot(132, sharex=axe, sharey=axe)
    axn = f.add_subplot(133, sharex=axe)
    #axj = f.add_subplot(144, sharex=axe)
    for index, i in enumerate(th.inj_index[ind]):
        lab=r't = {:.2f} s'.format(th.t[i])
        x=th.rho[i,:]
        axi.plot(x, th.nb_FKP_vars['pi'][i,:]*1e-3, col[index], lw=2, label=lab)
        axe.plot(x, th.nb_FKP_vars['pe'][i,:]*1e-3, col[index], lw=2, label=lab)
        axn.plot(x, th.nb_FKP_vars['n'][i, :]/th.kin_vars['ne'][i,:]*100., col[index], lw=2, label=lab)
        #axj.plot(x, th.nb_FKP_vars['nbcd'][i,:]*1e-3, col[index], lw=2, label=lab)
    au.limit_labels(axe,r'$\rho$', r'$P_e$ [$kW/m^3$]','')
    au.limit_labels(axi,r'$\rho$', r'$P_i$ [$kW/m^3$]','')
    au.limit_labels(axn,r'$\rho$', r'$n_f/n_e$ [%]','')
    #au.limit_labels(axj,r'$\rho$', r'j shielded [$kA/m^2$]','')

    axn.legend(bbox_to_anchor=(1.05, 0.65), loc=2)
    #f.text(0.01, 0.01, th.fname)

    f.tight_layout()
    f.subplots_adjust(right=0.8)
    plt.show()   


def plot_NBCD(th, time=[0]):
    """Plot current drive

    Makes a plot with NBCD (shielded and not), the shielding factor and the efficiency

    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        None    
    Note:        
    """
    au.common_style()

    plt.rc('legend', fontsize=15)

    try:
        th.nb_FKP_vars['nbcd'].mean()
    except:
        print("No nbcd data")
        return

    ind = tu._time_to_ind(th.t, time, inj_index=th.inj_index)
    f = plt.figure(figsize=(10,8))
    axj  = f.add_subplot(224)
    axcd = f.add_subplot(221)
    axsh = f.add_subplot(222)
    axeff = f.add_subplot(223)

    for ax in [axcd, axsh, axeff]:
        if time[0]!=0:
            for i, el in enumerate(ind):
                ax.axvline(x=th.t[th.inj_index[el]], color=col[i], lw=1.8, linestyle='--')

    axcd.plot(th.t[th.inj_index], th.nbcd[th.inj_index]*1e-3, 'k', lw=2.3, label='SH.')
    axcd.plot(th.t[th.inj_index], th.unbcd[th.inj_index]*1e-3, 'k:', lw=2.3, label='UNSH.')   
    axsh.plot(th.t[th.inj_index], 1.-th.shield[th.inj_index], 'k', lw=2.3)
    axeff.plot(th.t[th.inj_index], th.eff[th.inj_index]*100., 'k', lw=2.3)
    for i, el in enumerate(ind):
        axj.plot(th.rho[el,:], th.nb_FKP_vars['nbcd'][el,:]*1e-3, \
                 col[i],lw=2.3,label=r't = {:.2f} s'.format(th.t[el]))
    au.limit_labels(axj, r'$\rho$', r'$j^{SH}$ [$kA/m^2$]')
    au.limit_labels(axcd, r't [s]', r'$I_{CD}$ [$kA$]')
    axcd.legend(loc='best')
    au.limit_labels(axsh, r't [s]', r'$1-I_{SH}/I_{UN}$')
    au.limit_labels(axeff, r't [s]', r' $\eta \left[\frac{10^{18} A}{W m^2}\right]$', r'NBCD efficiency')
    axeff.set_ylim([0,1.])
    axj.legend(bbox_to_anchor=(-2.35,1.65), loc=2)
    #f.text(0.01, 0.01, th.runid)
    f.tight_layout();
    f.subplots_adjust(left=0.22, right=0.9)

    plt.show()


def plot_powerbalance(th, time, ls):
    """ plot pow balance

    Plots the power balance with cumulative plot
   
    Parameters:
        time (arr) : array with the times   where to plot the lines
    Attributes:
        None    
    Note:
    """   
    au.common_style()

    col, col2, styles, my_cmap, dpi = au.define_colors()

    #################################
    # defining vars to plot
    #################################
    ind = tu._time_to_ind(th.t, time, inj_index=th.inj_index)
    x,y,yH,yfrac = th._init_pbalance_vars()

    labels = ['Abs.', 'CX', 'O. L.', 'Shine-thr.']
    xlabel = r'Time [s]'
    ylabel = r'Power [MW]'
    
    #################################
    # Plot
    #################################
    f=plt.figure(figsize=(15,5))
    axfrac = f.add_subplot(122)
    axabs = f.add_subplot(121)
    au._cumulative_plot(x,y,labels, xlabel, ylabel, col2, ax=axabs, title='')

    if len(time)!=1:
        for i, el in enumerate(ind):
            axfrac.axvline(x=th.t[th.inj_index[el]], color='k', lw=2.3, linestyle='--')       
            axabs.axvline(x=th.t[th.inj_index[el]], color='k', lw=2.3, linestyle='--')      
    f.subplots_adjust(right=0.8)        

    for i, el in enumerate(yfrac):
        axfrac.plot(x, el*100., color=col2[i], lw=2.3, label=labels[i], linestyle=ls)
        
    au.limit_labels(axfrac, r'Time [s]', r'Fraction [%]')        

    axfrac.legend(loc='best')

    f.tight_layout()
    
