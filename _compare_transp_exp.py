#!/usr/bin/env python


import pytransp.classes.transp_exp as te
import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.integrate as integrate
import glob, sys

def _conv_exp_psi2phi(t, time):
    """Converts psi 2 phi

    Converts in between coordinates using the following relation (phi=toroidal flux, psi=poloidal flux)
    psi = int(1/q, phi)
    
    phi=int(q,psi)
    
    Paramters:
    None
    Arguments:
    None
    """
    q=t.getNode(r'\results::q_psi')
    _tq=q.getDimensionAt(1).data()
    _ind=np.argmin(_tq-time<0.)
    _psiq=q.getDimensionAt(0).data()
    _psiq=_psiq/float(np.max(_psiq))
    q=q.data()[_ind,:]
    locpsi = _psiq
    phi = integrate.cumtrapz(q,locpsi)
    phi = np.concatenate([[0], phi])
    phi = phi/max(phi)
    phip=interp.interp1d(locpsi, phi, fill_value=0.)
    # it returns the parameter of toroidal flux in poloidal flux
    return phip

if len(sys.argv) >= 3:
    shot = sys.argv[1]
    time = float(sys.argv[2])
    if len(sys.argv) == 4:
        tid = sys.argv[3]
else:
    print("Please give as input shot number, time, id (OPTIONAL)")
    print("Files will be looked for in '/home/vallar/tr_client/TCV/<<SHOT>>")
    print('\n e.g. \n compare_transp_exp.py 62124 1.3 (V01) \n')
    sys.exit()

folder='/home/vallar/tr_client/TCV/'+str(shot)+'/'
try:
    tid
    fname=folder+str(shot)+str(tid)+'.CDF'
except NameError:
    files = [f for f in glob.glob(folder + "*.CDF")]
    _np=0
    for i,f in enumerate(files):
        _n = f[-6:-4]
        if int(_n)<_np:
            continue
        else:
            _np=int(_n)
            ind_f=i
    try:
        ind_f
    except NameError:
        print('Files do not exist in '+folder)
        sys.exit()
        
    fname=files[ind_f]        
print('Getting CDF file '+fname)

t = mds.Tree('tcv_shot', int(shot))
output = te.transp_exp(fname)
data = {'EXP':{}, 'TRANSP':{}}
phip=_conv_exp_psi2phi(t,time)
#####################
# EXPERIMENTAL DATA
#####################

#vloop
vl_s = t.getNode('\magnetics::vloop')
vl_t = vl_s.getDimensionAt(0).data()
vl = vl_s.data()[0,:]
data['EXP']['vl'] = np.array([vl_t, vl])
#wmhd
wmhd_s = t.getNode(r'\results::total_energy:foo')
wmhd = wmhd_s.data()
wmhd_t = wmhd_s.getDimensionAt(0).data()
data['EXP']['wmhd'] = np.array([wmhd_t, wmhd])


#==============================================================
#beta
betat_s = t.getNode(r'\results::beta_tor')
betat = betat_s.data()
betat_t = betat_s.getDimensionAt(0).data()
#normalized beta = beta*a*btor/ip
# '100*beta/ip*1e6*b0*a_minor
rmax_s = t.getNode(r'\results::r_max_psi')
rmax = rmax_s.data()[:,-1]
rmin_s = t.getNode(r'\results::r_min_psi')
rmin = rmin_s.data()[:,-1]
print('Got stuff for betan')

beta_s = t.getNode(r'\results::beta_tor')
beta = beta_s.data()
betat_t = beta_s.getDimensionAt(0).data()
data['EXP']['betaN'] = np.array([betat_t, beta])
#==============================================================



ibs_s = t.getNode(r'\tcv_shot::top.results.ibs:ibs')
ibs_t = ibs_s.getDimensionAt(0).data()
ibs   = -1.*ibs_s.data()
data['EXP']['ibs'] = np.array([ibs_t, -1.*ibs])

pe_s=t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:pe')
rho_s=t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:rho')
pe=pe_s.data()
rho=rho_s.data()
pe_t = pe_s.getDimensionAt(0).data()
ind=np.argmin(pe_t-time<0.)
pe = pe[:,ind]

ne_s=t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:ne')
ne=ne_s.data()
ne_t = ne_s.getDimensionAt(0).data()
ind=np.argmin(ne_t-time<0.)
ne = ne[:,ind]

te_s=t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:te')
te=te_s.data()
te_t = te_s.getDimensionAt(0).data()
ind=np.argmin(te_t-time<0.)
te = te[:,ind]

ttexp=te_t[ind]

we_s = t.getNode(r'\tcv_shot::top.results.thomson.profiles.auto:we')
we   = te_s.data()
we_t = te_s.getDimensionAt(0).data()

ti_s=t.getNode(r'\tcv_shot::top.results.cxrs.proffit:ti')
ti=ti_s.data()
rho_ti = ti_s.getDimensionAt(0).data()
ti_t = ti_s.getDimensionAt(1).data()
ind=np.argmin(ti_t-time<0.)
ti = ti[ind,:]
rho_ti = rho_ti[ind,:]
_ind=rho_ti<1.
ti=ti[_ind]; rho_ti=rho_ti[_ind]

####################
# Converting rho of thomson (poloidal norm flux) to rho of TRANSP (toroidal norm flux
####################
tflux=phip(rho**2)
rhot=tflux**0.5

data['EXP']['pe'] = np.array([rhot, pe])
data['EXP']['ne'] = np.array([rhot, ne])
data['EXP']['te'] = np.array([rhot, te])
#data['EXP']['we'] = np.array([we_t, we])

tflux=phip(rho_ti**2)
rhot_ti=tflux**0.5
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
data['TRANSP']['n0'] = np.array([tt, output.n0_tot[:,-1]*1e6])
#data['TRANSP']['n0'] = np.array([[0],[0]])
ind=np.argmin(tt-time<0.)
ttransp=tt[ind]
data['TRANSP']['ne'] = np.array([output.rho[ind, :], output.kin_vars['ne'][ind, :]])
data['TRANSP']['te'] = np.array([output.rho[ind, :], output.kin_vars['te'][ind, :]])
data['TRANSP']['we'] = np.array([tt, output.we])

data['TRANSP']['pe'] = np.array([output.rho[ind, :], data['TRANSP']['ne'][1,:]*data['TRANSP']['te'][1,:]*1.602e-19])

#np.array([output.rho[ind, :], output.file.variables['PMHDT_IN'][ind, :]])

data['TRANSP']['ti'] = np.array([output.rho[ind, :], output.kin_vars['ti'][ind, :]])
#data['TRANSP']['ni'] = np.array([output.rho[ind, :], output.kin_vars['ni'][ind, :]])


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
            
    if indk==0:
        ax.legend(loc='best')
    if max(data['TRANSP'][key][0,:])<1 and min(data['TRANSP'][key][0,:])>0 or max(data['TRANSP'][key][0,:])>1:
        ax.set_xlabel(r'$\rho_{TOR}$')
    else:
        ax.set_xlabel(r't (s)')
    ax.grid('on')
f.suptitle(shot+'  t = '+str(time))
f.tight_layout()
plt.show()
