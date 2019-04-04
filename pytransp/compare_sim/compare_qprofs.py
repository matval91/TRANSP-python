import transp_output as to
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import define_colors
colors,_,_,_,_ = define_colors()

def compare_qprofs(fnames=['/home/vallar/TCV/58823/58823V68.CDF','/home/vallar/TCV/58823/58823V69.CDF'],\
                labels=['',''], time=1.25):

    f=plt.figure()
    ax=f.add_subplot(111)
    if labels[0]=='':
        labels = [fname[-12:-4] for fname in fnames]
    for i, fname in enumerate(fnames):
        try:
            o=to.transp_output(fname)
        except:
            print('Error! No file named '+fname)
            continue
        ind = np.argmin(o.t-time<0.)
        ax.plot(o.rho[ind,:], o.file.variables['Q'][ind,:], label=labels[i])

    ax.legend(loc='best')
    f.tight_layout()
    plt.show()
