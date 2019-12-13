#!/usr/bin/env python

import pytransp.classes.transp_heating as th
import pytransp.plot.plot_nb as plot_nb

import sys
   
def main(fnames=['../examples/65052V01.cdf'], time=1.):
    """
    """
    ths= [th.transp_heating(i) for i in fnames]; 
    print(ths[0].fname)
    plot_nb.plot_deposition_multishot(ths, time)

if len(sys.argv) >=3:
    time = float(sys.argv[1])
    fnames = sys.argv[2:]
    main(fnames, time)
else:
    print("Please give as input the time of interest and a/some CDF shot files")
    print("")
    print('\n e.g. \n plot_nubeam.py 1. 65052V01.CDF \n')
    sys.exit()
