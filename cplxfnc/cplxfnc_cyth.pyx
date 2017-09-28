cimport cython

cdef extern from "stdbool.h":
    ctypedef char bool

cdef extern from "../cplxfnc_clib/cplxfnc.hpp" namespace "cplxfnc":
    int zeta(double complex s, double complex a, double complex *res, double tol,
             unsigned int limit, bool verbose);


def py_zeta(double complex s, double complex a, double tol=1e-16, unsigned int limit=5, bool verbose=False):
    cdef double complex res
    cdef int status
    status = zeta(s, a, &res, tol, limit, verbose)
    if status != 0:
        raise RuntimeError("Zeta error")
    return res