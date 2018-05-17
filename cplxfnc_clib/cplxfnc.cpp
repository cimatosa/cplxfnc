/*
 *  MIT License
 *
 *  Copyright (c) 2017 Richard Hartmann
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */

#include "cplxfnc.hpp"

#include "acb.h"
#include "arb.h"
#include "arf.h"
#include "acb_hypgeom.h"
#include "flint/flint.h"

static_assert(__FLINT_VERSION == 2, "requires flint version 2");

#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>

extern "C" {
  void libcplxfnc_is_present(void) {}
}

namespace cplxfnc {

// ##################################################
// ##     Hurwitz Zeta function
// ##################################################

std::complex<double> zeta(std::complex<double> s, std::complex<double> a)
{
    return zeta(s, a, 1e-16, 5, false);
}

std::complex<double> zeta(std::complex<double> s, std::complex<double> a, double tol, unsigned int limit, bool verbose, unsigned int init_prec)
{
    std::complex<double> res;
    int status = zeta(s, a, &res, tol, limit, verbose, init_prec);
    if (status) {
        std::ostringstream oss;
        oss << "LIMIT ERROR: zeta s=" << s << " and a=" << a;
        throw std::runtime_error(oss.str());
    } else {
        return res;
    }
}

int zeta(std::complex<double> s, std::complex<double> a, std::complex<double> * res, double tol,
         unsigned int limit, bool verbose, unsigned int init_prec)
{
    acb_t _z, _s, _a;
    acb_init(_z); acb_init(_s); acb_init(_a);
    acb_set_d_d(_s, s.real(), s.imag());
    acb_set_d_d(_a, a.real(), a.imag());

    unsigned int prec = init_prec;
    unsigned int c = 1;
    double zeta_re, zeta_im;
    slong err_bits, err_bits_ref;
    err_bits_ref = slong(log2(tol));

    while (1) {
        acb_hurwitz_zeta(_z, _s, _a, prec);

        err_bits =  acb_rel_error_bits(_z);
        if (verbose) {
            std::cerr << std::setprecision(1) << std::fixed <<
            "zeta(s, a) with s=" << s << " and a=" << a << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
        }        
        
        if (err_bits <= err_bits_ref) {
            zeta_re = arf_get_d(arb_midref(acb_realref(_z)), ARF_RND_NEAR);
            zeta_im = arf_get_d(arb_midref(acb_imagref(_z)), ARF_RND_NEAR);
            acb_clear(_z); acb_clear(_s); acb_clear(_a);
            *res = std::complex<double>(zeta_re, zeta_im);
            return 0;
        }
        prec *= 2;
        c += 1;
        if (c > limit) {
            std::cerr << "\nERROR: limit (" << limit << ") reached\n" <<
            std::setprecision(1) << std::fixed <<
            "zeta(s, a) with s=" << s << " and a=" << a << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
            return -1;
        }
    }
}

// ##################################################
// ##     incomplete upper gamma function
// ##################################################

std::complex<double> gamma_inc(std::complex<double> s, std::complex<double> z)
{
    return gamma_inc(s, z, 1e-16, 5, false);
}

std::complex<double> gamma_inc(std::complex<double> s, std::complex<double> z, double tol,
                               unsigned int limit, bool verbose, unsigned int init_prec)
{
    std::complex<double> res;
    int status = gamma_inc(s, z, &res, tol, limit, verbose, init_prec);
    if (status) {
        std::ostringstream oss;
        if (status == -1) {
            oss << "LIMIT ERROR: gamma_inc s=" << s << " and z=" << z;
            throw std::runtime_error(oss.str());
        } else if (status == -2) {
            oss << "VALUE ERROR: gamma_inc, if Re(s) < 0 then z must not be zero!";
            throw std::runtime_error(oss.str());
        } else {
            oss << "gamma_inc unknown error: error code: " << status;
            throw std::runtime_error(oss.str());
        }
    } else {
        return res;        
    }
}

int gamma_inc(std::complex<double> s, std::complex<double> z, std::complex<double> * res, double tol,
         unsigned int limit, bool verbose, unsigned int init_prec)
{
    if ((s.real() < 0) && (z.real() == 0) && (z.imag() == 0)){
        if (verbose) {
            std::cerr << "ERROR: inc gamma value error!\n" <<
            "if Re(s) < 0 then z must not be zero!\n";
        }
        return -2;
    }
    acb_t _z, _s, _res;
    acb_init(_z); acb_init(_s); acb_init(_res);
    acb_set_d_d(_s, s.real(), s.imag());
    acb_set_d_d(_z, z.real(), z.imag());

    unsigned int prec = init_prec;
    unsigned int c = 1;
    double res_re, res_im;
    slong err_bits, err_bits_ref;
    err_bits_ref = slong(log2(tol));

    while (1) {
        acb_hypgeom_gamma_upper(_res, _s, _z, 0, prec);   //
        
        err_bits =  acb_rel_error_bits(_res);
        if (verbose) {
            std::cout << std::setprecision(1) << std::fixed <<
            "gamma(s, z) with s=" << s << " and z=" << z << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
        }
        
        if (err_bits <= err_bits_ref) {
            res_re     = arf_get_d(arb_midref(acb_realref(_res)), ARF_RND_NEAR);
            res_im     = arf_get_d(arb_midref(acb_imagref(_res)), ARF_RND_NEAR);
            acb_clear(_res); acb_clear(_s); acb_clear(_z);
            *res = std::complex<double>(res_re, res_im);
            return 0;
        }
        c += 1;
        if (c > limit) {
            std::cerr << "\nERROR: limit (" << limit << ") reached\n" <<
            std::setprecision(1) << std::fixed <<
            "gamma(s, z) with s=" << s << " and z=" << z << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
            return -1;
        }
        prec *= 2;

    }
}

// ##################################################
// ##     Asymptotic series for the confluent hypergeometric function
// ##     see: http://arblib.org/hypergeometric.html#algorithms-hypergeometric-asymptotic-confluent
// ##################################################


int u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z, std::complex<double> * res, double tol,
            unsigned int limit, bool verbose, unsigned int init_prec)
{
    acb_t _a, _b, _z, _res;
    acb_init(_a); acb_init(_b); acb_init(_z); acb_init(_res);
    acb_set_d_d(_a, a.real(), a.imag());
    acb_set_d_d(_b, b.real(), b.imag());
    acb_set_d_d(_z, z.real(), z.imag());
    
    unsigned int prec = init_prec;
    unsigned int c = 1;
    double res_re, res_im;
    slong err_bits, err_bits_ref;
    err_bits_ref = slong(log2(tol));
    
    if (not acb_hypgeom_u_use_asymp(_z, -err_bits_ref)) {
        std::cerr << "ERROR: u_asymp can not be evaluated for the given tolerence, this is a property of u_asymp!\n" <<
        "z:" << z << " tol:" << tol << " err_bits_ref:" << err_bits_ref << "\n" <<
        "acb_hypgeom_u_use_asymp(z, -err_bits_ref) failed\n"
        "increase z or decrease tol!\n";
        return -2;
    }
    
    while (1) {
        acb_hypgeom_u_asymp(_res, _a, _b, _z, -1, prec);   // n=-1 -> choose n automatically 
        
        err_bits =  acb_rel_error_bits(_res);
        if (verbose) {
            std::cout << std::setprecision(1) << std::fixed <<
            "u_asymp(a, b, z) with a=" << a << " and b=" << b << " and z=" << z << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
        }
        
        if (err_bits <= err_bits_ref) {
            res_re     = arf_get_d(arb_midref(acb_realref(_res)), ARF_RND_NEAR);
            res_im     = arf_get_d(arb_midref(acb_imagref(_res)), ARF_RND_NEAR);
            acb_clear(_res); acb_clear(_a); acb_clear(_b); acb_clear(_z);
            *res = std::complex<double>(res_re, res_im);
            return 0;
        }
        c += 1;
        if (c > limit) {
            std::cerr << "\nERROR: limit (" << limit << ") reached\n" <<
            std::setprecision(1) << std::fixed <<
            "u_asymp(a, b, z) with a=" << a << " and b=" << b << " and z=" << z << std::endl <<
            "internal prec: " << prec << std::endl <<
            "tol (bits)     : " << err_bits_ref << std::endl <<
            "rel_err (bits) : " << err_bits << std::endl;
            return -1;
        }
        prec *= 2;

    }
}

std::complex<double> u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z, double tol,
                             unsigned int limit, bool verbose, unsigned int init_prec)
{
    std::complex<double> res;
    int status = u_asymp(a, b, z, &res, tol, limit, verbose, init_prec);
    if (status) {
        std::ostringstream oss;
        if (status == -1) {
            oss << "LIMIT ERROR: u_asymp(a, b, z) with a=" << a << " and b=" << b << " and z=" << z;
            throw std::runtime_error(oss.str());
        } else if (status == -2) {
            oss << "VALUE ERROR: u_asymp can not be evaluated for the given (z, tolerence). This is a property of u_asymp! Increase z or decrease tol!";
            throw std::runtime_error(oss.str());
        } else {
            oss << "u_asymp unknown error: error code: " << status;
            throw std::runtime_error(oss.str());
        }
    } else {
        return res;
    }
}

std::complex<double> u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z)
{
    return u_asymp(a, b, z, 1e-16, 5, false);
}

} /* namespace cplxfnc */
