cimport cython

cdef extern from "stdbool.h":
    ctypedef char bool

cdef extern from "../cplxfnc_clib/cplxfnc.hpp" namespace "cplxfnc":
    double complex zeta(double complex s, double complex a, double tol, unsigned int limit, bool verbose) except +
    double complex gamma_inc(double complex s, double complex z, double tol, unsigned int limit, bool verbose) except +
    double complex u_asymp(double complex a, double complex b, double complex z, double tol, unsigned int limit, bool verbose) except +

def py_zeta(double complex s, double complex a, double tol=1e-16, unsigned int limit=5, bool verbose=False):
    return zeta(s, a, tol, limit, verbose)

def py_gamma_inc(double complex s, double complex z, double tol=1e-16, unsigned int limit=5, bool verbose=False):
    return gamma_inc(s, z, tol, limit, verbose)
    
def py_u_asymp(double complex a, double complex b, double complex z, double tol=1e-16, unsigned int limit=5, bool verbose=False):
    return u_asymp(a, b, z, tol, limit, verbose)    
