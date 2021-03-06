import transp_fbm
import matplotlib.pyplot as plt
from utils.plot_utils import limit_labels, common_style
common_style()

#file1 = input('Insert file 1: ')
#file2 = input('Insert file 2: ')

file1='/home/vallar/TCV/58823/58823V69_fi_1.cdf'
file2='/home/vallar/TCV/58823/58823V73_fi_2.cdf'

o1 = transp_fbm.fbm(file1)
o2 = transp_fbm.fbm(file2)

# f = plt.figure(figsize=(8,8))
# axne = f.add_subplot(211)
# axTe = f.add_subplot(212, sharex=axne)
# o1.plot_input_1d(axne=axne, axTe=axTe)
# axne.legend(loc='best')
# axTe.legend(loc='best')
# o2.plot_input_1d(axne=axne, axTe=axTe, ls='--')
# limit_labels(axne, r'Time [s]', r'$\langle n \rangle$ [$1/m^3$]','' )
# limit_labels(axTe, r'Time [s]', r'$\langle T \rangle$ [$eV$]','' )
# f.tight_layout()

# f=plt.figure(figsize=(10,8))
# axfrac = f.add_subplot(111)
# x,y,yH,yfrac = o1._init_pbalance_vars()
# o1._plot_pbalance_frac(x, yfrac, axfrac, ylim=[0,70])
# x,y,yH,yfrac = o2._init_pbalance_vars()
# o2._plot_pbalance_frac(x, yfrac, axfrac, ylim=[0,70], ls='--', leg=0)
# f.text(0.15, 0.95,'Shot #'+str(shot1)+' (on-axis,  solid )', fontsize=16, color='k')
# f.text(0.6, 0.95,'Shot #'+str(shot2)+' (off-axis, dotted)', fontsize=16, color='k')
# f.tight_layout()
# plt.subplots_adjust(top=0.9)

f = plt.figure()
axE = f.add_subplot(111)
#axp = f.add_subplot(212, sharex=axp)
axE.plot(o1.dict_dim['E']*1e-3, o1.f_spacep_int/o1.norm, 'k-', lw=2.3, label=r'No EC')
axE.plot(o2.dict_dim['E']*1e-3, o2.f_spacep_int/o2.norm, 'r-', lw=2.3, label=r'Full EC')
axE.set_xlim(0, 25)
limit_labels(axE,r'Energy [keV]',r'f [1/keV]','')
f.tight_layout()
axE.legend(loc='best')

o1._integrate_spaceE()
o2._integrate_spaceE()

f = plt.figure()
axE = f.add_subplot(111)
#axp = f.add_subplot(212, sharex=axp)
axE.plot(o1.dict_dim['pitch'], o1.f_spaceE_int, 'k-', lw=2.3, label=r'No EC')
axE.plot(o2.dict_dim['pitch'], o2.f_spaceE_int, 'r-', lw=2.3, label=r'Full EC')
axE.set_xlim(-1., 1.)
limit_labels(axE,r'$\xi=v_\parallel/v$',r'f ','')
f.tight_layout()
axE.legend(loc='best')

plt.show()
