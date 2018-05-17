import mpmath as mp
import cplxfnc as cf
import numpy as np

import os
import sys
from contextlib import contextmanager

@contextmanager
def stdouterr_redirected(to=os.devnull):
    '''
    taken from
    https://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python

    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd_out = sys.stdout.fileno()
    fd_err = sys.stderr.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1


    def _redirect_stdout(to):
        sys.stdout.close()
        os.dup2(to.fileno(), fd_out)        # fd_out writes to 'to' file
        sys.stdout = os.fdopen(fd_out, 'w') # Python writes to fd_out

    def _redirect_stderr(to):
        sys.stderr.close()
        os.dup2(to.fileno(), fd_err)        # fd_err writes to 'to' file
        sys.stderr = os.fdopen(fd_err, 'w') # Python writes to fd_err

    with os.fdopen(os.dup(fd_out), 'w') as old_stdout:
        with os.fdopen(os.dup(fd_err), 'w') as old_stderr:
            with open(to, 'w') as file:
                _redirect_stdout(to=file)
                _redirect_stderr(to=file)
            try:
                yield # allow code to be run with the redirected stdout
            finally:
                _redirect_stdout(to=old_stdout) # restore stdout. buffering and flags such as
                _redirect_stderr(to=old_stderr)


def cplx_rand(re_low, re_high, im_low, im_high):
    return np.random.rand() * (re_high-re_low) + re_low + 1j*(np.random.rand() * (im_high-im_low) + im_low)

def test_zeta(n=50, tol=1e-16):
    np.random.seed(0)
    mp.mp.dps = 64
    for i in range(n):
        s = cplx_rand(-5, 5, -5, 5)
        a = cplx_rand(-5, 5, -5, 5)
        z = cf.zeta(s, a, tol=tol)
        z_mp = mp.zeta(s, a)
        assert (abs(z - complex(z_mp)))/abs(complex(z_mp)) < tol

def test_gamma_inc(n=50, tol=1e-16):
    mp.mp.dps = 64
    for i in range(n):
        s = cplx_rand(-5, 5, -5, 5)
        z = cplx_rand(-5, 5, -5, 5)
        g = cf.gamma_inc(s, z, tol)
        g_mp = mp.gammainc(s, z)
        assert (abs(g - complex(g_mp))) / abs(complex(g_mp)) < tol


def test_uasymp():
    """
        here we test the equivalence of Gamma(-s, z) exp(z) z^(s+1) = u_asymp(s+1,s+1,z)
        and compare with accurate results from mpmath

        for 0 < s < 10 and z = +- 100

        which allows to use

        Gamma(-s, z) exp(z) z^(s+1)    for abs(z) < 100

        and

        u_asymp(s+1,s+1,z)    for abs(z) >= 100
    """
    mp.mp.dps = 64

    eps = sys.float_info.epsilon

    def u_asymp_mp_real(s, z):
        s = mp.mpf(s)
        z = mp.mpf(z)
        U = mp.gammainc(z=-s, a=z) * mp.exp(z) * z ** (1 + s)
        return float(U.real)

    def gamma_inc_mp(s, z):
        s = mp.mpf(s)
        z = mp.mpf(z)
        U = mp.gammainc(z=s, a=z)
        return complex(U)

    for s in np.linspace(0, 10, 75):

        # the case z = w/wc < 0
        abs_z = 100
        pow_err_fac = 250

        u = cf.u_asymp(s+1, s+1, abs_z)
        u_ref = u_asymp_mp_real(s, abs_z)
        d = abs(u - u_ref)
        assert d < 1e-15, "U_ASYMP: u:{}, u_ref:{}, d:{}, s:{}, abs_z:{}".format(u, u_ref, d, s, abs_z)

        ig = cf.gamma_inc(-s, abs_z)
        ig_mp = gamma_inc_mp(-s, abs_z)

        d = abs(ig - ig_mp)/abs(ig_mp)
        assert d < eps, "gamma_inc seems inaccurate"

        ez = np.exp(abs_z)
        ez_mp = float(mp.exp(abs_z))
        d = abs(ez - ez_mp) / abs(ez_mp)
        assert d < eps, "np.exp seems inaccurate"

        p = abs_z**(s+1)
        p_mp = float(mp.mpf(abs_z)**(mp.mpf(s)+1))
        d = abs(p - p_mp)

        d_max = pow_err_fac*(s+1)*abs_z**s*eps
        # power_err_fac account of the inaccuracy of the pow implementation of the c math library
        assert d < d_max, "pow seems inaccurate {}, {}".format(d, d_max)

        u = ig * np.exp(abs_z) * abs_z**(s+1)
        d_max = abs(ig * np.exp(abs_z)) * d_max
        d = abs(u - u_ref)
        assert d < d_max, "IG: u:{}, u_ref:{}, d:{}, s:{}, z:{}".format(u, u_ref, d, s, abs_z)

        # z = w/wc > 0
        z = abs_z

        u = cf.u_asymp(s + 1, s + 1, -z)
        u_ref = u_asymp_mp_real(s, -z)

        d = abs(u - u_ref)
        assert d < 1e-15, "U_ASYMP: u:{}, u_ref:{}, d:{}, s:{}, abs_z:{}".format(u, u_ref, d, s, z)

        ig_eps = cf.gamma_inc(-s, -(z + 1j * eps))
        ig = cf.gamma_inc(-s, -z)
        if ig_eps.imag * ig.imag < 0:
            ig = ig.conjugate()

        u = ig * np.exp(-z) * z ** (s + 1) * np.exp(-1j*np.pi*(s+1))

        d = abs(u - u_ref)
        assert d < pow_err_fac*eps, "IG: u:{}, u_ref:{}, d:{}, s:{}, z:{}".format(u, u_ref, d, s, abs_z)


    try:
        with stdouterr_redirected():
            u = cf.u_asymp(2, 2, -15)
    except RuntimeError:
        pass
    else:
        assert False, "expected RuntimeError"

    try:
        with stdouterr_redirected():
            u = cf.u_asymp(2, 2, 30)
    except RuntimeError:
        pass
    else:
        assert False, "expected RuntimeError"


if __name__ == "__main__":
    test_zeta(10)
    test_gamma_inc(10)
    test_uasymp()