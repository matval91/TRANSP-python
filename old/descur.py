import ctypes as ct
import numpy as np
import socket
class DESCUR:


    def descur_fit(self, r_surf, z_surf, mom_order, tok='AUG'):

#2D fit Total LSQR? fit 
        if (socket.gethostname()[:4]) == 'lac5':
            libd = '/home/vianello/pythonlib/tcv/nbh/descur/lib/descur_idl.so'
        else:
            libd = '/home/vianello/pythonlib/tcv/nbh/descurlac/lib/descur_idl.so'
            
        libdsc = ct.cdll.LoadLibrary(libd)

        momtype = 4
        nphi = 1
        c_nmom  = ct.c_long(mom_order)
        c_ok    = ct.c_long(0)
        c_nphi  = ct.c_long(nphi)
        c_niter = ct.c_long(4000)
        c_nstep = ct.c_long(100)
        c_nfp   = ct.c_long(1)
        c_ftol  = ct.c_double(1.E-9)
        c_pexp  = ct.c_double(4.E0)
        c_qexp  = ct.c_double(4.E0)
        nphi1   = nphi+1
        nphi3   = nphi+3
        c_rmnaxis = (ct.c_float * nphi1)()
        c_zmnaxis = (ct.c_float * nphi1)()
        c_rbc = (ct.c_double * mom_order * nphi3)()
        c_zbc = (ct.c_double * mom_order * nphi3)()
        c_rbs = (ct.c_double * mom_order * nphi3)()
        c_zbs = (ct.c_double * mom_order * nphi3)()

        moments = np.zeros((mom_order, momtype))
        nrz = len(r_surf)
        c_nrz = ct.c_long(nrz)
        if (nrz <=  0):
            return moments
        c_rin = (ct.c_double * nrz)()
        c_zin = (ct.c_double * nrz)()
        for j_nrz in range(nrz):
            c_rin[j_nrz] = r_surf[j_nrz]
            c_zin[j_nrz] = z_surf[j_nrz]

        _nmom    = ct.byref(c_nmom)
        _nrz     = ct.byref(c_nrz)
        _nphi    = ct.byref(c_nphi)
        _niter   = ct.byref(c_niter)
        _nstep   = ct.byref(c_nstep)
        _nfp     = ct.byref(c_nfp)
        _ftol    = ct.byref(c_ftol)
        _pexp    = ct.byref(c_pexp)
        _qexp    = ct.byref(c_qexp)
        _rin     = ct.byref(c_rin)
        _zin     = ct.byref(c_zin)
        _ok      = ct.byref(c_ok)
        _rbc     = ct.byref(c_rbc)
        _rbs     = ct.byref(c_rbs)
        _zbc     = ct.byref(c_zbc)
        _zbs     = ct.byref(c_zbs)
        _rmnaxis = ct.byref(c_rmnaxis)
        _zmnaxis = ct.byref(c_zmnaxis)
        if (min(z_surf) < -1.0) and (tok=='AUG'):
            print('Unacceptable field lines')
        else:
            status = libdsc.curve3d_(_nmom, _nrz, _nphi, _niter, _nstep, _nfp, \
                                     _ftol, _pexp, _qexp, _rin, _zin, _ok, \
                                     _rbc, _zbs, _rbs, _zbc, _rmnaxis, _zmnaxis)

            moments[:, 0] = c_rbc[1]
            moments[:, 1] = c_rbs[1]
            moments[:, 2] = c_zbc[1]
            moments[:, 3] = c_zbs[1]

        return moments


    def descur_fit_fast(self, r_surf, z_surf, mom_order):

# LSQR decoposition of the R,Z flux surface coordinates to Fourier moments
# r_surf - (nt,nrho,nangle) or (nrho,nangle) or (nangle)
# higher mom_order is neccessary! 

        nangle = np.size(r_surf, -1)

        rb = np.fft.rfft(r_surf)[..., :mom_order]/nangle
        zb = np.fft.rfft(z_surf)[..., :mom_order]/nangle

        moments = np.zeros(r_surf.shape[:-1] + (mom_order, 4))

        moments[...,  :mom_order, 0] = np.real(rb)
        moments[...,  :mom_order, 1] = np.imag(rb)
        moments[...,  :mom_order, 2] = np.real(zb)
        moments[...,  :mom_order, 3] = np.imag(zb)
        moments[..., 1:mom_order, :] *= 2

        return moments
