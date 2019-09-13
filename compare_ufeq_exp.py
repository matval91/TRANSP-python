import ufiles_omf as uf
import numpy as np
import matplotlib.pyplot as plt
import MDSplus as mds


shot=63234
conn = mds.Connection('tcvdata.epfl.ch')
conn.openTree('tcv_shot', shot)
liuqe_flavour='LIUQE.M'
_d={'label':'', 'x':[], 'y':[]}
signals={'Fvac':dict.copy(_d),
         'L2B':dict.copy(_d),
         'torflux':dict.copy(_d),
         'polflux':dict.copy(_d),
         'psi_axis':dict.copy(_d),
         'F':dict.copy(_d)}
signals['Fvac']['label']='rbtor_vac'
signals['L2B']['label']='lambda'
signals['torflux']['label']='tor_flux_tot'
signals['polflux']['label']='psi_surf'
signals['psi_axis']['label']='psi_axis'
signals['F']['label']='rbtor_rho'

for i in signals.keys():
    _s = conn.get(r"tcv_eq('"+signals[i]['label']+"', '"+liuqe_flavour+"')")
    if i=='F':
        _time=conn.get('dim_of(tcv_eq("'+signals[i]['label']+'", "'+liuqe_flavour+'"),1)').data()
        ind=np.where(np.logical_and(_time>1.29, _time<1.31))[0]
    if i!='torflux':
        s = _s.data()
        s_t=conn.get('dim_of(tcv_eq("'+signals[i]['label']+'", "'+liuqe_flavour+'"),0)').data()
    else:
        s = _s.data()[:,-1]
        s_t=conn.get('dim_of(tcv_eq("'+signals[i]['label']+'", "'+liuqe_flavour+'"),1)').data()
    
    signals[i]['x'] = s_t
    signals[i]['y'] = s

#signals['F']['y'] *= 0.996 #correction
    

f=uf.RU('/home/vallar/63234_transp/OMF63234.RBZ')

fig=plt.figure(); ax=fig.add_subplot(111)
ax.plot(f.values['X0'], f.fvalues, 'kx', label='UF')
ax.plot(signals['Fvac']['x'], signals['Fvac']['y'], 'r', label='exp')
ax.set_xlabel(r't [s]'); ax.set_ylabel(r'R*btor (vac) [T*m]')
plt.legend(loc='best')
plt.tight_layout()
print('Rbtor exp/uf',np.mean(signals['F']['y'])/np.mean(f.fvalues))

# f=uf.RU('/home/vallar/63234_transp/OMF63234.L2B')

# fig=plt.figure(); ax=fig.add_subplot(111)
# ax.plot(f.values['X0'], f.fvalues, 'kx', label='UF')
# ax.plot(signals['L2B']['x'], signals['L2B']['y'], 'r', label='exp')
# ax.set_xlabel(r't [s]'); ax.set_ylabel(r'Lambda')
# plt.legend(loc='best')

# plt.show()

f=uf.RU('/home/vallar/63234_transp/OMF63234.TRF')
f=uf.RU('/home/vallar/63234_input/AA63234.TRF')

fig=plt.figure(); ax=fig.add_subplot(111)
ax.plot(f.values['X0'], -1.*f.fvalues, 'k', label='UF')
ax.plot(signals['torflux']['x'], signals['torflux']['y'], 'r', label='exp')
ax.set_xlabel(r't [s]'); ax.set_ylabel(r'torflux')
plt.legend(loc='best')
plt.tight_layout()

print('torflux exp/uf',np.mean(signals['torflux']['y'])/np.mean(f.fvalues))

plt.show()
plt.grid('on')


f=uf.RU('/home/vallar/63234_transp/OMF63234.PLF')
#f=uf.RU('/home/vallar/63234_input/AA63234.PLF')

fig=plt.figure(); ax=fig.add_subplot(111)
ax.plot(f.values['X0'], f.fvalues, 'k', label='UF')
ax.plot(signals['polflux']['x'], signals['polflux']['y']/(2*np.pi), 'r', label='exp')
ax.set_xlabel(r't [s]'); ax.set_ylabel(r'polflux')
plt.legend(loc='best')
print('polflux exp/uf',np.mean(signals['polflux']['y']/(2.*np.pi))/np.mean(f.fvalues))
plt.tight_layout()

f=uf.RU('/home/vallar/63234_transp/OMF63234.GRB')
f=uf.RU('/home/vallar/63234_input/AA63234.GRB')

fig=plt.figure(); ax=fig.add_subplot(111)
ax.plot(f.values['X0'], f.fvalues*-1., 'k', label='UF')
ax.plot(signals['F']['x'], signals['F']['y'][ind,:].T, 'r', label='exp')
ax.set_xlabel(r'rho'); ax.set_ylabel(r'R Bt [T*m]')
#plt.legend(loc='best')
#print('polflux exp/uf',np.mean(signals['polflux']['y']/(2.*np.pi))/np.mean(f.fvalues))
plt.tight_layout()

plt.show()
