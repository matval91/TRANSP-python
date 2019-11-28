
import sys
import pytransp.classes.transp_deposition as td
def main(fname= '../examples/65052V01_fi_1.cdf', fname_surf=''):
    """
    """
    dep = td.absorption(fname, fname_surf)
    # plotting
    dep.plot_XYpart()
    dep.plot_RZ()
    dep.plot_Epitch()
    dep.plot_rhopitch()
    return

if len(sys.argv) ==2:
    fname = sys.argv[1]
    r,z,phi,E,xi = main(fname)
elif len(sys.argv) ==3:
    fname = sys.argv[1]
    fname_surf = sys.argv[2]
    main(fname, fname_surf)
else:
    print("Please give as input a birth profile")
    print('\n e.g. \n nubeam_plot_deposition_profile.py ../examples/65052V02_birth.cdf1 \n')
    sys.exit()
    
    
    
    

