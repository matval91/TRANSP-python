#!/usr/bin/env python

import pytransp.compare_wexp.scale_taup as st
import sys
import numpy as np
import matlab.engine

if len(sys.argv)==1:
    print("Please give as input the shots id and requested time. The simulations will be searched in /home/vallar/tr_client/TCV/<<shot>>")
    print('\n e.g. \n compare_transp_sim.py 62124V01 62124V02 1.3 \n')
    sys.exit()
