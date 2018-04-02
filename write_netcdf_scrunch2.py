import numpy as np
import os
import netCDF4 as nc

def write_netcdf_scrunch2(geo_struct, sig_struct):
        """
        Writing file movie.cdf, as suggested by M. Gorelenkova on 3/30/18, to give to scrunch for reconstruction.
        Just needs the toroidal flux, not the moments.

        allocate(xary(inx),zary(inz),psirz(inx,inz,intimes),psi_efit(inx))
        allocate(times(intimes),gzero(intimes),xplas(intimes),apl(intimes))
        allocate(gary(inpsi,intimes),qprof2(inpsi,intimes))
        allocate(ggprime(inpsi,intimes))
        allocate(pres(inpsi,intimes),presp(inpsi,intimes))
        allocate(psirad(inpsi,intimes))
        allocate(rbnd(intheta,intimes),zbnd(intheta,intimes))
        allocate(inpsis(intimes),inbnds(intimes),iminx(inz),imaxx(inz))
        allocate(geq_filnams(intimes))

        ! read data from "movie.cdf"
        
        ! these are time-dependant
        call cdf_read(icdf,'times',times)
        call cdf_read(icdf,'gzero',gzero)
        call cdf_read(icdf,'xplas',xplas)
        call cdf_read(icdf,'apl',apl)
        ! these are (psi, times)
        call cdf_read(icdf,'xsv2',psirad)
        call cdf_read(icdf,'gary',gary)
        call cdf_read(icdf,'qprof2',qprof2)
        call cdf_read(icdf,'ggprime',ggprime)
        call cdf_read(icdf,'pres',pres)
        call cdf_read(icdf,'presp',presp)
        ! these are (nR, nz)
        call cdf_read(icdf,'xary',xary)
        call cdf_read(icdf,'zary',zary)
        call cdf_read(icdf,'psi',psirz)
        ! these are (theta, times)
        call cdf_read(icdf,'rbnd',rbnd)
        call cdf_read(icdf,'zbnd',zbnd)
        """
        if os.path.isfile('movie.cdf'):
            print("Removing old movie.cdf")
            os.remove('movie.cdf')

        
        geo_struct._calc_RBtor();  geo_struct._calc_qprof()
        geo_struct._get_pressure(); geo_struct._get_torflux(); geo_struct.mag_axis()
        torFlux = geo_struct.eq.rz2phinorm(geo_struct.R, geo_struct.Z, geo_struct.time_u, make_grid=True)
        torFlux[np.isnan(torFlux)] = 1.
        torFlux = np.transpose(torFlux, (1,2,0)) #put it in shape (nR, nZ, time)

        sig_struct._read_1d()
        
        
        f = nc.Dataset('movie.cdf','w', format='NETCDF4_CLASSIC')

        # time variables
        t = f.createDimension('t', geo_struct.nt)
        dummy = np.zeros(geo_struct.nt)
        times = f.createVariable('times', np.float, ('t',))
        gzero = f.createVariable('gzero', np.float, ('t',))
        xplas = f.createVariable('xplas', np.float, ('t',))
        apl   = f.createVariable('apl',   np.float, ('t',))
        times[:] = geo_struct.time_u
        gzero[:] = geo_struct.RBT[:, 0]
        xplas[:] = geo_struct.r_axis
        apl[:] = sig_struct._univec['IP']['spline'](sig_struct.time_u)

        # psi variables (dummy)
        psi = f.createDimension('psi', 41)
        dummyarr = np.full((41, geo_struct.nt), 1., dtype=float)
        te2     = f.createVariable('te2'    , np.float, ('psi', 'psi',))
        xsv2    = f.createVariable('xsv2'   , np.float, ('t', 'psi',))
        gary    = f.createVariable('gary'   , np.float, ('t', 'psi',))
        qprof2  = f.createVariable('qprof2' , np.float, ('t', 'psi',))
        ggprime = f.createVariable('ggprime', np.float, ('t', 'psi',))
        pres    = f.createVariable('pres'   , np.float, ('t', 'psi',))
        presp   = f.createVariable('presp'  , np.float, ('t', 'psi',))
        xsv2[:] = geo_struct.phi
        gary[:]    = geo_struct.RBT
        qprof2[:]  = geo_struct.qprof
        ggprime[:] = geo_struct.FFP
        pres[:] = geo_struct.pressure
        presp[:] = geo_struct.pprime
        
        # nR, nZ
        R = f.createDimension('R', len(geo_struct.R))
        Z = f.createDimension('Z', len(geo_struct.Z))
        xary = f.createVariable('xary', np.float, ('R',))
        zary = f.createVariable('zary', np.float, ('Z',))
        psi  = f.createVariable('psi' , np.float, ('t','R','Z',))
        xary[:] = geo_struct.R; zary[:]=geo_struct.Z
        psi[:] = torFlux

        # theta, times
        theta = f.createDimension('theta', geo_struct.n_the)
        rlcfs = f.createVariable('rbnd', np.float, ('t', 'theta',))
        zlcfs = f.createVariable('zbnd', np.float, ('t', 'theta',))
        geo_struct._eqflux_LCFS()
        rlcfs[:] = np.transpose(geo_struct.r_plot_LCFS)
        zlcfs[:] = np.transpose(geo_struct.z_plot_LCFS)

        #close cdf file
        f.close()
        
