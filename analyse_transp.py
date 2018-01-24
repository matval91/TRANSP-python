import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import patches

colors = ['k','r','b','g','m']
colors = ['m','g','b','r','k','c']
def cumulative_plot(x,y,labels, xlabel, ylabel):
    f  = plt.figure()
    ax = f.add_subplot(111)
    tmpy=np.zeros(len(x))
    for i, el in enumerate(y):
        if i!=len(y)-1:
            tmpy +=el
            ax.plot(x, tmpy, colors[i], lw=3., label=labels[i])
        else:
            ax.plot(x, el, colors[i], lw=3., label=labels[i])

        if i==0:
            y0 = np.zeros(len(x))
            ax.fill_between(x, y0, tmpy, color=colors[i])
        elif i!=len(y)-1:
            ax.fill_between(x, tmpy-el, tmpy, color=colors[i])
        elif i==len(y)-1:
            ax.fill_between(x, tmpy, el, color=colors[i])
            
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax2=ax.twinx()
    ticks1 = ax.get_yticks()
    ax2.set_yticks(ticks1)
    ax2.set_yticklabels(ticks1/np.max(el))
    ax2.set_ylim(ax.get_ylim())
    ax2.set_ylabel('Fraction (%)')
    f.legend(loc='upper left')

def compute_Ec(A, Aj, Zj, ne, Te, nj):
    """
    Ec = 14.8*Te*(A^{3/2}/ne*sum(nj*Zj/Aj))^2/3
    """
    sum_facts = [nj[i]*Zj[i]**2/Aj[i] for i in range(len(Aj))]
    summ = np.sum(np.asarray(sum_facts), axis=0)
    second_factor = (A**1.5/ne*summ)**0.66
    Ec_prof = 14.8*Te*second_factor
    
    return Ec_prof

def compute_Gi(E0, Ec):
    """
    Gi = Ec/E0 int[0, E0/Ec]dy/(1+y**1.5)
    """
    gi = np.zeros(np.shape(Ec), dtype=float)
    for time_ind in range(np.shape(Ec)[0]):
        timeslice = Ec[time_ind,:]
        if E0[time_ind]==0:
            gi[time_ind,:] = np.zeros(np.shape(Ec)[1])
        else:
            for rho_ind in range(np.shape(Ec)[1]):
                xarray = np.linspace(0,E0[time_ind]/timeslice[rho_ind], num=100, dtype=float)
                gi[time_ind,rho_ind] = timeslice[rho_ind]/E0[time_ind]*np.trapz(1/(1+xarray**1.5), xarray)
    return gi

def compute_tau(A, Te, ne, E0,Ec):
    """
    # A plasma
    tau_s=ts/3*ln(1+(E/EC)**1.5)
    ts=6.27e8*A Te[eV]**1.5/(Z**2 ne(cm^-3) lnLambda)
    valid for Z=1
    """
    if np.mean(ne)/1e18 >1:
        ne = ne*1e-6
    else:
        ne=ne
    ts=6.27e8*A[0]*np.divide(np.power(Te,1.5),(1.*17*ne))
    E0=np.transpose(np.tile(E0, (np.shape(Ec)[1], 1)))
    taus=ts/3.*(1+(E0/Ec)**1.5)
    return taus


directory ='/home/vallar/MST1/T16/57720_transp/'
directory = '/home/vallar/tr_client/TCV/58823/B01/'
#fname = '57720A01.CDF'
fname = '58823B01.CDF'
infile = directory+fname

fff = nc.Dataset(infile)

darea = fff.variables['DAREA'][:]
dvol  = fff.variables['DVOL'][:]

rho = fff.variables['X'][:] # rho with time dependence
time = fff.variables['TIME'][:]

Abeam = 2
E0 = fff.variables['EINJAV_D'][:]
Aplasma = [2,12]
Zplasma = [1,6]
# """
# Variation of kinetic variables i.e. density and temperature (and beam energy) with beam
# """
# # time-rho variables
ne = fff.variables['NE'][:]
ni = fff.variables['NI'][:]
nimp = fff.variables['NIMP'][:]
# nfast = fff.variables['NB01_TOT']

te = fff.variables['TE'][:]
ti = fff.variables['TI']
# ebeam = fff.variables['EINJ01_E1']
# f=plt.figure()
# ax1=f.add_subplot(111)
# ax1.set_title(r'Shot #57720')
# ax1.plot(time, ne[:,0], 'k', label=r'n_e')
# ax1.plot(time, ni[:,0], 'k--', label=r'n_i')
# ax1.plot(time, nfast[:,0], 'k:', label=r'n_{fast}')
# ax1.set_xlabel(r'Time (s)')
# ax1.set_ylabel(r'Density ($cm^{-3}$)')
# ax2 = ax1.twinx()
# ax2.plot(time, te[:,0]*1e-3, 'r', label=r'$t_{e}$')
# ax2.plot(time, ti[:,0]*1e-3, 'r--', label=r'$t_{i}$')
# ax2.plot(time, ebeam[:]*0.1*1e-3, 'r:', label=r'$0.1 \times E_{fast}$')
# ax2.set_ylabel(r'Temperature (keV)')
# ax2.tick_params('y', colors='r')
# f.legend(loc='best')

# """
# Beam power absorption
# """
# time-rho variables
pinj = fff.variables['PBINJ_D']
#pshine = fff.variables['PBSHINE_D']
#orbloss = fff.variables['BPLIM_D']
##pcx = fff.variables['BPCX0_D'][:]  # D BEAM CX SCE POWER (EXT)
##pcx = fff.variables['BPCI0_D'][:] # D BEAM CX SCE POWER (INT)
#pcx  = fff.variables['BPCXX_D'][:] # D BEAM POWER TO CX (EXT)
#pcx += fff.variables['BPCXI_D'][:] # D BEAM POWER TO CX (INT)
##pcx += fff.variables['BPCRX_D'][:] # FAST ION CX RECAPTURE (EXT)
##pcx += fff.variables['BPCRI_D'][:] # FAST ION CX RECAPTURE (INT)
pion = fff.variables['BPTI'][:] 
pel = fff.variables['BPTE'][:]
#ptherm = fff.variables['BPTH_D']
#sump = pshine[:]+orbloss[:]+pcx[:]+pion+pel
#
#ynbi=[pshine[:]*1e-6,orbloss[:]*1e-6, pcx[:]*1e-6,(pion[:]+pel[:])*1e-6,\
#      ptherm[:]*1e-6, pinj[:]*1e-6]
#labnbi=[r'P_{ST}', r'P_{OL}', r'P_{CX}', r'P_{ABS}',r'P_{th}', r'P_{inj}']
#cumulative_plot(time, ynbi, labnbi, 'Time (s)', 'Power (MW)')
#
#f=plt.figure()
#ax1=f.add_subplot(111)
#ax1.set_title(r'Shot #57720')
#ax1.plot(time, pinj[:]*1e-6, 'k', label=r'P_{inj}')
#ax1.plot(time, pshine[:]*1e-6, 'k--', label=r'P_{ST}')
#ax1.plot(time, orbloss[:]*1e-6, 'k-.', label=r'P_{OL}')
#ax1.plot(time, pcx[:]*1e-6, 'b', label=r'P_{CX}')
#ax1.plot(time, (pion[:]+pel[:])*1e-6, 'b', label=r'P_{ABS}')
#ax1.plot(time, ptherm[:]*1e-6, 'k:', label = r'P_{THERM}')
#ax1.plot(time, sump*1e-6, 'r-')
#
#
#ax1.set_xlabel(r'Time (s)')
#ax1.set_ylabel(r'Power (MW)')
#f.legend(loc='best')
# """
# Beam power absorption profiles
# """
#pbe_prof = fff.variables['PBE'][:]
#pbi_prof = fff.variables['PBI'][:]
#power_e_prof = [np.multiply(dvol[t,:], pbe_prof[t,:]) for t in range(len(time))]
#power_i_prof = [np.multiply(dvol[t,:], pbi_prof[t,:]) for t in range(len(time))]
#Ec = compute_Ec(Abeam, Aplasma, Zplasma, ne, te, [ni,nimp])
#Gi = compute_Gi(E0,Ec)
#taus = compute_tau(Aplasma, te, ne, E0, Ec)
#
#f = plt.figure()
#ax1 = f.add_subplot(421)
#ax2 = f.add_subplot(422, sharex=ax1) #121 plots on first row, first column, 122 plots on first row second column
#ax3 = f.add_subplot(423, sharex=ax1)
#ax4 = f.add_subplot(424, sharex=ax1)
#ax5 = f.add_subplot(425, sharex=ax1)
#ax6 = f.add_subplot(426, sharex=ax1)
#ax7 = f.add_subplot(427, sharex=ax1)
#
#for i in np.linspace(80, 150, num=3, dtype=int):
#    ax1.plot(rho[i,:], pbe_prof[i,:]*1e3, lw=2.3, label=pel[i])
#    ax2.plot(rho[i,:], pbi_prof[i,:]*1e3, lw=2.3, label=pion[i])
#    ax3.plot(rho[i,:], power_e_prof[i]*1e-3, lw=2.3, label=time[i])
#    ax4.plot(rho[i,:], power_i_prof[i]*1e-3, lw=2.3, label=time[i])
#    ax5.plot(rho[i,:], Ec[i,:]*1e-3, lw=2.3, label=time[i])
#    ax5.plot(rho[i,:], np.tile(E0[i]*1e-3, len(rho[i,:])), 'k',lw=2.5, label = 'E0')
#
#    ax6.plot(rho[i,:], Gi[i,:], lw=2.3, label=time[i])
#    ax6.plot(rho[i,:], 1-Gi[i,:], lw=2.3,linestyle='--', label=time[i])
#    ax7.plot(rho[i,:], taus[i,:]*1e3, lw=2.3,linestyle='--', label=time[i])
#    
#ax1.legend(loc='best'); ax2.legend(loc='best')
#
#ax1.set_xlabel(r'$\rho$'); ax2.set_xlabel(r'$\rho$'); ax3.set_xlabel(r'$\rho$');
#ax4.set_xlabel(r'$\rho$'); ax5.set_xlabel(r'$\rho$'); ax6.set_xlabel(r'$\rho$');
#ax7.set_xlabel(r'$\rho$')
#ax1.set_ylabel(r'Power density (kW/m^3)'); ax2.set_ylabel(r'Power density (kW/m^3)'); ax3.set_ylabel(r'Power (kW)');
#ax4.set_ylabel(r'Power (kW)'); ax5.set_ylabel(r'$E_c$ (keV)'); ax6.set_ylabel(r'$G_{i}$');
#ax7.set_ylabel(r'$\tau_s $ (ms)')
#
#ax1.set_title(r'Power density to electrons');ax2.set_title(r'Power density to ions'); ax3.set_title(r'Power to electrons'); 
#ax4.set_title(r'Power to ions'); ax5.set_title(r'Critical Energy'); ax6.set_title(r'Fraction to ions')
#f.tight_layout()
# """
# Current
# """
jtot_prof =  fff.variables['CUR'][:]
joh_prof = fff.variables['CUROH'][:]
nbcd_prof = fff.variables['CURB'][:]
jboot_prof = fff.variables['CURBS'][:]
jtot = np.zeros(len(time), dtype=float)
nbcd = np.zeros(len(time), dtype=float)
jboot = np.zeros(len(time), dtype=float)
joh = np.zeros(len(time), dtype=float)

for i in range(len(time)):
    nbcd[i] = np.dot(nbcd_prof[i,:], darea[i,:])
    jboot[i] = np.dot(jboot_prof[i,:], darea[i,:])
    joh[i] = np.dot(joh_prof[i,:], darea[i,:])
    jtot[i] = np.dot(jtot_prof[i,:], darea[i,:])

# j_exp = nbcd+joh+jboot
# #li = fff.variables['LI_1'][:]
f=plt.figure()
ax1=f.add_subplot(111)
ax1.set_title(r'Shot #57720')
ax1.plot(time, jtot[:]*1e-3, 'r', label=r'J_{TOT}', linewidth=3)
ax1.plot(time, joh[:]*1e-3, 'k', label=r'J_{OH}')
ax1.plot(time, nbcd[:]*1e-3, 'k-.'  , label=r'J_{CD}')
ax1.plot(time, jboot[:]*1e-3, 'k--', label=r'J_{BS}')

ax1.set_xlabel(r'Time (s)')
ax1.set_ylabel(r'Current (kA)')
f.legend(loc='best')


f=plt.figure()
ax1=f.add_subplot(111)
ax1.set_title(r'Shot #57720')
for i in np.linspace(80, 150, num=3, dtype=int):
    ax1.plot(rho[i,:], nbcd_prof[i,:], label=time[i])
    
ax1.set_ylabel('NBCD density (kA/m$^2$)')
 #ax2 = ax1.twinx()
 #ax2.plot(time, li, 'r', label=r'l_i')
 #ax2.set_ylabel(r'Internal inductance')
 #ax2.tick_params('y', colors='r')

f.legend(loc='best')

# """
# Beta, internal inductance variation
# """
# 1D: BETAE      units:
# 1D: BETAR      units:
# 1D: BETAI      units:
# 1D: BETAT      units:
# <type 'netCDF4._netCDF4.Variable'>
# float32 BETAE(TIME)
#     units:
#     long_name: ELECTRON BETA (POLOIDAL)
# unlimited dimensions:
# current shape = (201,)
# filling off

# <type 'netCDF4._netCDF4.Variable'>
# float32 BETAR(TIME)
#     units:
#     long_name: ROTATION BETA (POLOIDAL)
# unlimited dimensions:
# current shape = (201,)
# filling off

# <type 'netCDF4._netCDF4.Variable'>
# float32 BETAI(TIME)
#     units:
#     long_name: THERMAL ION BETA POLOIDAL
# unlimited dimensions:
# current shape = (201,)
# filling off

# <type 'netCDF4._netCDF4.Variable'>
# float32 BETAT(TIME)
#     units:
#     long_name: TOTAL BETA(POLOIDAL)
# unlimited dimensions:
# current shape = (201,)
# filling off

# <type 'netCDF4._netCDF4.Variable'>
# float32 LIO2(TIME)
#     units:
#     long_name: INDUCTANCE (LI/2)
# unlimited dimensions:
# current shape = (201,)
# filling off
li = fff.variables['LIO2'][:]
betae = fff.variables['BETAE'][:]
betai = fff.variables['BETAI'][:]
betat = fff.variables['BETAT'][:]
pinjnorm = pinj/np.max(pinj)
pec = fff.variables['PECIN'][:]
pecnorm = pec/np.max(pec)
tinj = time[pinjnorm>0.1]
tec = time[pecnorm>0.1]

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(time, li, lw=2.3, label=r'$l_i/2$')
ax.plot(time, betae, lw=2.3, label=r'$\beta_e$')
ax.plot(time, betai, lw=2.3, label=r'$\beta_i$')
ax.plot(time, betat, lw=2.3, label=r'$\beta_t$')
#ax.plot(time, pinjnorm, lw=2.3, label=r'$P/P_{MAX}$')
pinj=patches.Rectangle(
        (np.min(tinj), 0.), np.max(tinj)-np.min(tinj), np.max(ax.get_ylim()),
        alpha=0.2, facecolor='red' )
pec=patches.Rectangle(
        (np.min(tec), 0.), np.max(tec)-np.min(tec), np.max(ax.get_ylim()),
        alpha=0.2, facecolor='green' )
ax.add_patch(pinj); ax.add_patch(pec)
f.legend(loc='best')
plt.show()
