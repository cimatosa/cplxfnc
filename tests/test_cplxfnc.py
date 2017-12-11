import mpmath as mp
import cplxfnc as cf
import numpy as np

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
        g = cf.zeta(s, z, tol)
        g_mp = mp.zeta(s, z)
        assert (abs(g - complex(g_mp))) / abs(complex(g_mp)) < tol


if __name__ == "__main__":
    test_zeta(10)
    test_gamma_inc(10)