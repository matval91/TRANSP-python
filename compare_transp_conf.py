#!/usr/bin/env python

import pytransp.classes.transp_exp as te
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.integrate as integrate
import glob, sys

if len(sys.argv) >= 3:
    shot = sys.argv[1]
    time = float(sys.argv[2])
    if len(sys.argv) == 4:
        tid = sys.argv[3]
else:
    print("Please give as input shot number, time, id (OPTIONAL)")
    print("Files will be looked for in '/home/vallar/NoTivoli/tr_client/TCV/<<SHOT>>")
    print('\n e.g. \n compare_transp_exp.py 62124 1.3 (V01) \n')
    sys.exit()

folder='/home/vallar//NoTivoli/tr_client/TCV/'+str(shot)+'/'
try:
    tid
    fname=folder+str(shot)+str(tid)+'.CDF'
except NameError:
    files = [f for f in glob.glob(folder + "*.CDF")]
    _np=0; ind_f=0;
    for i,f in enumerate(files):
        _n = f[-6:-4]
        if int(_n)<_np:
            continue
        else:
            _np=int(_n)
        ind_f=i
    fname=files[ind_f]

print('Getting CDF file '+fname)

#t = mds.Tree('tcv_shot', int(shot))
conn = mds.Connection('tcvdata.epfl.ch')
conn.openTree('tcv_shot', int(shot))

output = te.transp_exp(fname)
data = {'EXP':{}, 'TRANSP':{}}
#####################
# EXPERIMENTAL DATA
#####################

#vloop
sig=r'\magnetics::vloop'
vl_s = conn.get(sig)
vl_t = conn.get('dim_of('+sig+',0)').data()
#vl_t = vl_s.getDimensionAt(0).data()
vl = vl_s.data()[0,:]
data['EXP']['vl'] = np.array([vl_t, vl])
#wmhd
sig=r'\results::total_energy:foo'
wmhd_s = conn.get(sig)
wmhd = wmhd_s.data()
wmhd_t = conn.get('dim_of('+sig+',0)').data()
data['EXP']['wmhd'] = np.array([wmhd_t, wmhd])


#==============================================================
#beta
sig=r'\results::beta_tor'
betat_s = conn.get(sig)
beta = betat_s.data()
betat_t = conn.get('dim_of('+sig+',0)').data()
data['EXP']['betaN'] = np.array([betat_t, beta])
#==============================================================

sig=r'\tcv_shot::top.results.ibs:ibs'
ibs_s = conn.get(sig)
ibs_t = conn.get('dim_of('+sig+',0)').data()
ibs   = ibs_s.data()
data['EXP']['ibs'] = np.array([ibs_t, -1.*ibs])

sig=r'\tcv_shot::top.results.conf:pe'
pe_s = conn.get(sig)
pe_t = conn.get('dim_of('+sig+',1)').data()
pe=pe_s.data()

sig = r'\tcv_shot::top.results.conf:rhotor'
rho_s = conn.get(sig)
rho = rho_s.data()
ind = np.argmin(pe_t-time<0.)
pe = pe[ind,:]
rho = rho[ind,:]

sig=r'\tcv_shot::top.results.conf:ne'
ne_s=conn.get(sig)
ne=ne_s.data()
ne_t = conn.get('dim_of('+sig+',1)').data()
ind=np.argmin(ne_t-time<0.)
ne = ne[ind,:]

sig=r'\tcv_shot::top.results.conf:te'
te_s=conn.get(sig)
te=te_s.data()
te_t = conn.get('dim_of('+sig+',1)').data()
ind=np.argmin(te_t-time<0.)
te = te[ind,:]

ttexp=te_t[ind]
sig=r'\tcv_shot::top.results.conf:ti'
ti_s=conn.get(sig)
ti=ti_s.data()
ti_t = conn.get('dim_of('+sig+',1)').data()
ind=np.argmin(ti_t-time<0.)
ti = ti[ind,:]

sig=r'\tcv_shot::top.results.conf:z_eff'
zeff_s=conn.get(sig)
zeff=zeff_s.data()
zeff_t = conn.get('dim_of('+sig+',0)').data()
####################
# Converting rho of thomson (poloidal norm flux) to rho of TRANSP (toroidal norm flux
####################
#tflux=phip(rho**2)
rhot=rho
rhot_ti=rho
data['zeff'] = np.array([zeff_t, zeff])
data['EXP']['pe'] = np.array([rhot, pe])
data['EXP']['ne'] = np.array([rhot, ne])
data['EXP']['te'] = np.array([rhot, te])
#data['EXP']['we'] = np.array([we_t, we])

data['EXP']['ti'] = np.array([rhot_ti, ti])
#data['EXP']['ni'] = np.array([rhot_ti, ni])

#N0 from ascii file
try:
    n0=np.loadtxt('/home/vallar/edge_density/'+str(shot)+'n0.dat')
    data['EXP']['n0'] = np.array([n0[:,0], n0[:,1]])
except:
    print('No n0 data')
    data['EXP']['n0'] = np.array([[0],[0]])

#####################
# TRANSP DATA
#####################
#output._calculate_scalars()
tt=output.t
data['TRANSP']['vl'] = np.array([tt, output.file.variables['VSURC']])
output._calculate_wmhd()
data['TRANSP']['wth'] = np.array([tt, output.wth])
data['TRANSP']['wmhd'] = np.array([tt, output.wtot])
data['TRANSP']['betaN'] = np.array([tt, output.file.variables['BTDIA']])
data['TRANSP']['ibs'] = np.array([tt, output.jboot])
data['TRANSP']['n0'] = np.array([tt, output.n0_tot[:,-1]])
#data['TRANSP']['n0'] = np.array([[0],[0]])
ind=np.argmin(tt-time<0.)
ttransp=tt[ind]
rho=output.rho[ind,:]
data['TRANSP']['ne'] = np.array([rho, output.kin_vars['ne'][ind, :]])
data['TRANSP']['te'] = np.array([rho, output.kin_vars['te'][ind, :]])
data['TRANSP']['wth'] = np.array([tt, output.wth])

data['TRANSP']['pe'] = np.array([output.rho[ind, :], data['TRANSP']['ne'][1,:]*data['TRANSP']['te'][1,:]*1.602e-19])

#np.array([output.rho[ind, :], output.file.variables['PMHDT_IN'][ind, :]])

data['TRANSP']['ti'] = np.array([rho, output.kin_vars['ti'][ind, :]])
#data['TRANSP']['ni'] = np.array([output.rho[ind, :], output.kin_vars['ni'][ind, :]])
data['TRANSP']['zeff'] = np.array([tt, output.file.variables['ZEFFI0']])

plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=10, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=10)
ncol=3; nrow=3

colexp ='k'
coltr  = 'b'
col = {'EXP':'k', 'TRANSP':'b'}
lab = {'EXP':'EXP (t={:2.2f})'.format(ttexp), 'TRANSP':'TRANSP (t={:2.2f})'.format(ttransp)}

f,axs=plt.subplots(nrow, ncol, figsize=[4*ncol,3*nrow])

for indk, key in enumerate(data['EXP'].keys()):
    if key=='zeff':
        continue
    _col=indk% ncol
    _row=int(indk/ncol)
    ax=axs[_row,_col]
    for kkk in ['EXP', 'TRANSP']:
        ax.plot(data[kkk][key][0,:], data[kkk][key][1,:], color=col[kkk], label=lab[kkk])
        ax.set_ylabel(key)

    if key=='wmhd':
        offset=input('Need an offset? It will be ADDED to experimental wmhd \n')
        if offset!=0:
            ax.plot(data['EXP'][key][0,:], data['EXP'][key][1,:]+offset, color=col['EXP'], ls='--', label='EXP w offset')
            ax.legend(loc='best')

        ax.plot(data['TRANSP']['wth'][0,:],data['TRANSP']['wth'][1,:] , color=col['TRANSP'], ls='-.')
    elif key=='vl':
        ax2=ax.twinx()
        ax2.plot(data['zeff'][0,:], data['zeff'][1,:], color=col['EXP'],ls='--')
        ax2.plot(data['TRANSP']['zeff'][0,:], data['TRANSP']['zeff'][1,:], color=col['TRANSP'], ls='--')
        ax.set_ylim([-2, 5.]); ax2.set_ylim([0., 5])
        ax2.set_ylabel(r'Z_{eff} (dotted)')
            
            
    if indk==0:
        ax.legend(loc='best')
    if max(data['TRANSP'][key][0,:])<=1 and min(data['TRANSP'][key][0,:])>0:
        ax.set_xlabel(r'$\rho_{TOR}$')
    else:
        ax.set_xlabel(r't (s)')
    ax.grid('on')
f.suptitle(shot+'  t = '+str(time))
f.tight_layout()

output.plot_powerbalance()

plt.show()
