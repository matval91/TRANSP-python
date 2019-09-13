#!/usr/bin/env python

import ufiles_omf as uf
import numpy as np
import matplotlib.pyplot as plt
import MDSplus as mds


shot=63234

shot=63704
folder='/home/vallar/missions/1953/63704/TRANSP/UFILES/'

#shot=63717
#folder='/home/vallar/missions/1953/63717/TRANSP/UFILES/'

#conn = mds.Connection('tcvdata.epfl.ch')
#conn.openTree('tcv_shot', shot)
t = mds.Tree('tcv_shot', int(shot))

#liuqe_flavour='LIUQE.M'
_d={'label':'', 't':[], 'rho':[], 'y':[], 'ufile':'', 'ylab':''}
signals={'ne':dict.copy(_d),
         'te':dict.copy(_d),
         'ti':dict.copy(_d)}

signals['ne']['label']=r'\tcv_shot::top.results.conf:ne'
signals['te']['label']=r'\tcv_shot::top.results.conf:te'
signals['ti']['label']=r'\tcv_shot::top.results.conf:ti'

signals['ne']['ufile']='OMF'+str(shot)+'.NER'
signals['te']['ufile']='OMF'+str(shot)+'.TER'
signals['ti']['ufile']='OMF'+str(shot)+'.TI2'

signals['ne']['ylab']=r'$n_e$'
signals['te']['ylab']=r'$t_e$'
signals['ti']['ylab']=r'$t_i$'

signals['ne']['fac']=1e6
signals['te']['fac']=1e3
signals['ti']['fac']=1e3

rho_s=t.getNode(r'\tcv_shot::top.results.conf:rhotor')
rho = rho_s.data()

for i in signals.keys():
    _s=t.getNode(signals[i]['label'])
    signals[i]['y'] = _s.data()
    signals[i]['t'] = _s.getDimensionAt(1).data()
    signals[i]['rho'] = rho

#ind=np.argmin(signals['ne']['t']-time<0.)


for i in signals.keys():
    print('Open ufile '+folder+signals[i]['ufile'] )
    f=uf.RU(folder+signals[i]['ufile']) #x0 is rho, x1 is time
    # plot time, n0
    x=f.values['X1']; y=f.fvalues[0,:]
    fig=plt.figure(); ax=fig.add_subplot(111)
    ax.plot(x,y*signals[i]['fac'], 'kx', label='UF')
    ax.plot(signals[i]['t'], signals[i]['y'][:,0], 'r', label='exp')
    ax.set_xlabel(r't [s]'); ax.set_ylabel(signals[i]['ylab'])
    plt.legend(loc='best')
    plt.tight_layout()
plt.show()
# # f=uf.RU('/home/vallar/63234_transp/OMF63234.L2B')

# # fig=plt.figure(); ax=fig.add_subplot(111)
# # ax.plot(f.values['X0'], f.fvalues, 'kx', label='UF')
# # ax.plot(signals['L2B']['x'], signals['L2B']['y'], 'r', label='exp')
# # ax.set_xlabel(r't [s]'); ax.set_ylabel(r'Lambda')
# # plt.legend(loc='best')

# # plt.show()

# f=uf.RU('/home/vallar/63234_transp/OMF63234.TRF')

# fig=plt.figure(); ax=fig.add_subplot(111)
# ax.plot(f.values['X0'], -1.*f.fvalues, 'k', label='UF')
# ax.plot(signals['torflux']['x'], signals['torflux']['y'], 'r', label='exp')
# ax.set_xlabel(r't [s]'); ax.set_ylabel(r'torflux')
# plt.legend(loc='best')
# plt.tight_layout()

# print('torflux exp/uf',np.mean(signals['torflux']['y'])/np.mean(f.fvalues))

# plt.show()
# plt.grid('on')


# f=uf.RU('/home/vallar/63234_transp/OMF63234.PLF')
# #f=uf.RU('/home/vallar/63234_input/AA63234.PLF')

# fig=plt.figure(); ax=fig.add_subplot(111)
# ax.plot(f.values['X0'], f.fvalues, 'k', label='UF')
# ax.plot(signals['polflux']['x'], signals['polflux']['y']/(2*np.pi), 'r', label='exp')
# ax.set_xlabel(r't [s]'); ax.set_ylabel(r'polflux')
# plt.legend(loc='best')
# print('polflux exp/uf',np.mean(signals['polflux']['y']/(2.*np.pi))/np.mean(f.fvalues))
# plt.tight_layout()

# f=uf.RU('/home/vallar/63234_transp/OMF63234.GRB')
# f=uf.RU('/home/vallar/63234_input/AA63234.GRB')

# fig=plt.figure(); ax=fig.add_subplot(111)
# ax.plot(f.values['X0'], f.fvalues*-1., 'k', label='UF')
# ax.plot(signals['F']['x'], signals['F']['y'][ind,:].T, 'r', label='exp')
# ax.set_xlabel(r'rho'); ax.set_ylabel(r'R Bt [T*m]')
# #plt.legend(loc='best')
# #print('polflux exp/uf',np.mean(signals['polflux']['y']/(2.*np.pi))/np.mean(f.fvalues))
# plt.tight_layout()

# plt.show()
