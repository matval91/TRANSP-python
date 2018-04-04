import numpy as np
import os
import netCDF4 as nc


def write_netcdf_scrunch2(g):
        """
        Writing file movie.cdf, as suggested by M. Gorelenkova on 3/30/18, to give to scrunch for reconstruction.
        Needs everything since it is recomputing the equilibrium with better tolerances

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

        ###===========================
        # Getting data
        ###===========================
        
        g._calc_RBtor();  g._calc_qprof()
        g._get_pressure(); g._get_torflux(); g.mag_axis()
        g._eqflux_LCFS()
        R_wall = [0.6, 1.16]
        z_wall = [-0.77, 0.77]
        #R = np.linspace(np.min(g.r_plot_LCFS), np.max(g.r_plot_LCFS), 256)
        #Z = np.linspace(np.min(g.z_plot_LCFS), np.max(g.z_plot_LCFS), 256)
        tim = g.eq.getTimeBase()
        ind_t = np.argmin(tim-1.<0.)
        torFlux = g.eq.getFluxGrid()[ind_t,:,:]
        torFlux = np.expand_dims(torFlux, axis=0)
        torFlux = np.transpose(torFlux, (1,2,0)) #put it in shape (nR, nZ, time)
        R  = np.linspace(R_wall[0], R_wall[1], np.shape(torFlux)[1])
        Z  = np.linspace(z_wall[0], z_wall[1], np.shape(torFlux)[0])

        ###===========================
        # writing data to netcdf
        ###===========================        
        f = nc.Dataset('movie.cdf','w', format='NETCDF3_CLASSIC')
        # time variables
        t = f.createDimension('t', g.nt)
        dummy = np.zeros(g.nt)
        times = f.createVariable('times', np.float, ('t',))
        gzero = f.createVariable('gzero', np.float, ('t',))
        xplas = f.createVariable('xplas', np.float, ('t',))
        apl   = f.createVariable('apl',   np.float, ('t',))
        times[:] = g.time_u
        gzero[:] = g.RBT[:, 0]
        xplas[:] = g.r_axis*100.
        apl[:] = g.eq.getIpMeas()[ind_t]

        # psi variables 
        psi = f.createDimension('psi', 41)
        dummyarr = np.full((41, g.nt), 1., dtype=float)
        te2     = f.createVariable('te2'    , np.float, ('psi', 'psi',))
        xsv2    = f.createVariable('xsv2'   , np.float, ('t', 'psi',))
        gary    = f.createVariable('gary'   , np.float, ('t', 'psi',))
        qprof2  = f.createVariable('qprof2' , np.float, ('t', 'psi',))
        ggprime = f.createVariable('ggprime', np.float, ('t', 'psi',))
        pres    = f.createVariable('pres'   , np.float, ('t', 'psi',))
        presp   = f.createVariable('presp'  , np.float, ('t', 'psi',))
        xsv2[:] = g.phi
        gary[:]    = g.RBT*100.
        qprof2[:]  = g.qprof
        ggprime[:] = g.FFP*100.*100.
        pres[:] = g.pressure
        presp[:] = g.pprime
        
        # nR, nZ
        Rnc = f.createDimension('R', len(R))
        Znc = f.createDimension('Z', len(Z))
        xary = f.createVariable('xary', np.float, ('R',))
        zary = f.createVariable('zary', np.float, ('Z',))
        psi  = f.createVariable('psi' , np.float, ('t','Z','R',))
        xary[:] = R*100.; zary[:]=Z*100.
        psi[:] = torFlux

        # theta, times
        theta = f.createDimension('theta', g.n_the)
        rlcfs = f.createVariable('rbnd', np.float, ('t', 'theta',))
        zlcfs = f.createVariable('zbnd', np.float, ('t', 'theta',))
        rlcfs[:] = np.transpose(g.r_plot_LCFS*100.)
        zlcfs[:] = np.transpose(g.z_plot_LCFS*100.)

        #close cdf file
        f.close()
        
