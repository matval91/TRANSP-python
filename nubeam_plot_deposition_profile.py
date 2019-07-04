import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import math

shot = '57722'
run = 'A02'
fname = shot+run+'_nubeam_birth.cdf1'
home='/home/vallar/'
#dir = home+'/NoTivoli/OUTPUT_trial_nubeam/'
#dir=home
#dir = '/home/vallar/NUBEAM_src/nubeam/nubeam_comp_exec/TEST_TCV/test/'
#fname = '/home/vallar/52252D01_nubeam_birth.cdf1'
#fname = '/home/vallar/52252B03_nubeam_birth.cdf1'
try:
    var = nc.Dataset(home+'/'+fname,'r').variables
    title = home+''+fname
except:
    var = nc.Dataset(dir+fname,'r').variables
    title = dir+fname
    
r=var['bs_r_D_MCBEAM'][:]*0.01
z=var['bs_z_D_MCBEAM'][:]*0.01
phi=var['bs_zeta_D_MCBEAM'][:]/180*math.pi

Rmin = 0.6 # in meters
Rmax = 1.16 #in meters
zmin=-0.77
zmax=0.77
phiall = np.linspace(0,2*math.pi,100)
r_box = [Rmin, Rmin, Rmax, Rmax]
z_box = [zmin, zmax, zmax, zmin]

f1 = plt.figure()
ax_xy = f1.add_subplot(111)
ax_xy.plot(r*np.cos(phi), r*np.sin(phi), 'x')
ax_xy.plot(Rmin*np.cos(phiall), Rmin*np.sin(phiall), '.')
ax_xy.plot(Rmax*np.cos(phiall), Rmax*np.sin(phiall), '.')
ax_xy.set_xlabel(r'X (m)')
ax_xy.set_ylabel(r'Y (m)')

f2 = plt.figure()
ax_rz = f2.add_subplot(111)
hb=ax_rz.hexbin(r,z,gridsize=100, cmap='inferno')
lb = f2.colorbar(hb, ax=ax_rz)
#ax_rz.plot(r, z, 'x')
ax_rz.plot(r_box, z_box, 'k-')
ax_rz.set_xlabel(r'R (m)')
ax_rz.set_ylabel(r'Z (m)')


ax_xy.set_title(title)
ax_rz.set_title(title)
plt.show()
