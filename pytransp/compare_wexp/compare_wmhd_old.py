import transp_output
import numpy as np
import matplotlib.pyplot as plt
colors = ['r','b','g', 'c', 'm']
shot=58823;
#shot=58832
if shot==58823:
    fname1 = '/home/vallar/TCV/58823/58823V65.CDF' #TAUP=5ms
    fname2 = '/home/vallar/TCV/58823/58823V66.CDF' #TAUP=5ms
    #fname3 = '/home/vallar/TCV/58823/58823V62.CDF' #TAUP=5ms
    fname4 = '/home/vallar/TCV/58823/58823V68.CDF' #TAUP=5ms
    fname5 = '/home/vallar/TCV/58823/58823V99.CDF' #TAUP=5ms

    fnames = [fname1,  fname2, fname4, fname5]
    exp_fname_w = '/home/vallar/TCV/58823/58823wmhd.dat'
    exp_fname_v = '/home/vallar/TCV/58823/58823vl.dat'
elif shot==58832:
    #fname1 = '/home/vallar/TCV/58823/58823O03.CDF'
    #fname2 = '/home/vallar/TCV/58823/58823O05.CDF'
    fname1 = '/home/vallar/TCV/58832/58832V19.CDF' #TAUP=5ms
    fname2 = '/home/vallar/TCV/58832/58832V20.CDF' #TAUP=5ms
    fname3 = '/home/vallar/TCV/58832/58832V22.CDF' #TAUP=5ms
    fname4 = '/home/vallar/TCV/58832/58832V23.CDF' #TAUP=5ms
    fname5 = '/home/vallar/TCV/58832/58832V24.CDF' #TAUP=5ms

    fnames = [fname1, fname2, fname3, fname4, fname5]
    exp_fname_w = '/home/vallar/TCV/58832/58832wmhd.dat'
    exp_fname_v = '/home/vallar/TCV/58832/58832vl.dat'
elif shot==62124:
    #fname1 = '/home/vallar/TCV/58823/58823O03.CDF'
    #fname2 = '/home/vallar/TCV/58823/58823O05.CDF'
    fname1 = '/home/vallar/TCV/62124A01.CDF' #TAUP=5ms


    fnames = [fname1]#, fname2, fname3, fname4, fname5]
    exp_fname_w = '/home/vallar/TCV/62124wmhd.dat'
    exp_fname_v = '/home/vallar/TCV/62124vl.dat'
vldata = np.loadtxt(exp_fname_v)
wmhddata = np.loadtxt(exp_fname_w)

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(wmhddata[:,0], wmhddata[:,1]*1e-3, 'k--', label='EXP')
fv=plt.figure()
axv=fv.add_subplot(111)
axv.plot(vldata[:,0], vldata[:,1], 'k--', label='EXP')
fz=plt.figure(); axz=fz.add_subplot(111);

for i, fname in enumerate(fnames):
    try:
        o=transp_output.output_1d(fname)
    except:
        continue
    o._calculate_wmhd()
    ax.plot(o.t, o.wtot*1e-3, colors[i], label=fname[-12:-4])
    ax.plot(o.t, o.wth*1e-3, colors[i], linestyle='--', label=fname[-12:-4])
    #axv.plot(o.t,o.file.variables['VSURC'][:], colors[i], label=fname[-12:-4])
    #axv.plot(o.t,o.file.variables['VSUR'][:], colors[i], label=fname[-12:-4])
    axz.plot(o.t, o.file.variables['ZEFFM'][:], colors[i], label=fname[-12:-4])
    axz.plot(o.t, o.file.variables['ZEFFC'][:], colors[i], label=fname[-12:-4])

ax.set_xlabel(r'Time (s)'); ax.set_ylabel(r'W$_{MHD}$ [kJ]')
ax.legend(loc='best')
f.tight_layout()

#axv.set_xlabel(r'Time (s)'); axv.set_ylabel(r'V$_{LOOP}$')
#axv.legend(loc='best')
#fv.tight_layout()
#
#axz.set_xlabel(r'Time (s)'); axz.set_ylabel(r'Z$_{EFF}$')
#axz.legend(loc='best')
#fz.tight_layout()

ax.grid('on'); #axv.grid('on'); axz.grid('on')
plt.show()
