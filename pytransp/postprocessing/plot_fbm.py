import pytransp.classes.transp_fbm as tfbm
import sys
   
def main(fname='../examples/65052V01_fi_1.cdf'):
    """
    """
    fbm=tfbm.transp_fbm(fname)
    fbm.plot_space()
    fbm.plot_Epitch()
    fbm.plot_spaceE()
    fbm.plot_spacep()
    return fbm

if len(sys.argv) ==2:
    fname = sys.argv[1]
    main(fname)
else:
    print("Please give as input a FBM function")
    print("")
    print('\n e.g. \n plot_fbm 65052V01_fi_1.cdf \n')
    sys.exit()

