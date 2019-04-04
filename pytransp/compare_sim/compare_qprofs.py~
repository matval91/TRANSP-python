import transp_output
import numpy as np
import matplotlib.pyplot as plt
colors = ['r','b','g', 'c', 'm']
shot=58823;
#shot=58832
if shot==58823:
    fname1 = '/home/vallar/TCV/58823/58823V35.CDF' #TAUP=5ms
    #fname2 = '/home/vallar/TCV/58823/58823V60.CDF' #TAUP=5ms
    fname3 = '/home/vallar/TCV/58823/58823V61.CDF' #TAUP=5ms
    fname4 = '/home/vallar/TCV/58823/58823V50.CDF' #TAUP=5ms
    fname5 = '/home/vallar/TCV/58823/58823V58.CDF' #TAUP=5ms

    fnames = [fname1,  fname3]#, fname4, fname5]

f=plt.figure()
ax=f.add_subplot(111)
t=0.9
for i, fname in enumerate(fnames):
    try:
        o=transp_output.output_1d(fname)
    except:
        continue
    ind = np.argmin(o.t-t<0.)
    ax.plot(o.file.variables['X'][ind,:], o.file.variables['Q'][ind,:], label=fname[-12:-4])

ax.legend(loc='best')
f.tight_layout()
plt.show()
