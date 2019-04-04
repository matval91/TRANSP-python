import pytransp.classes.transp_output as to
import matplotlib.pyplot as plt
from utils.ascot_utils import limit_labels

def compare_sim(files=['/home/vallar/TCV/58823/58823V68.CDF','/home/vallar/TCV/58823/58823V99.CDF' ],\
                labels=['','']):
    """
    """
    data = {}
    for i in range(len(files)):
        data[str(i)] = {}
        data[str(i)]['fname'] = files[i]
        data[str(i)]['label'] = files[labels[i]]
        data[str(i)]['file']  = to.transp_output(files[i])

    #f = plt.figure(figsize=(8,8))
    #axne = f.add_subplot(211)
    #axTe = f.add_subplot(212, sharex=axne)
    o1.plot_input_prof()
    axne.legend(loc='best')
    axTe.legend(loc='best')
    f = plt.figure(figsize=(8,8))
    o2.plot_input_prof()
    #limit_labels(axne, r'Time [s]', r'$\langle n \rangle$ [$1/m^3$]','' )
    #limit_labels(axTe, r'Time [s]', r'$\langle T \rangle$ [$eV$]','' )
    f.tight_layout()

    f=plt.figure(figsize=(10,8))
    axfrac = f.add_subplot(111)
    x,y,yH,yfrac = o1._init_pbalance_vars()
    o1._plot_pbalance_frac(x, yfrac, axfrac, ylim=[0,70])
    x,y,yH,yfrac = o2._init_pbalance_vars()
    o2._plot_pbalance_frac(x, yfrac, axfrac, ylim=[0,70], ls='--', leg=0)
    f.text(0.15, 0.95,'Shot #'+str(shot1)+' (on-axis,  solid )', fontsize=16, color='k')
    f.text(0.6, 0.95,'Shot #'+str(shot2)+' (off-axis, dotted)', fontsize=16, color='k')
    f.tight_layout()
    plt.subplots_adjust(top=0.9)

    f = plt.figure(figsize=(8,8))
    axp = f.add_subplot(211)
    axf = f.add_subplot(212, sharex=axp)
    o1._plot_deposition_1d(axp=axp, axf=axf)
    axp.legend(loc='best')
    axf.legend(loc='best')
    o2._plot_deposition_1d(axp=axp, axf=axf)
    limit_labels(axp,r'Time [s]',r'P[kW]','')
    limit_labels(axf,r'Time [s]',r'%','')
    f.tight_layout()


    plt.show()
