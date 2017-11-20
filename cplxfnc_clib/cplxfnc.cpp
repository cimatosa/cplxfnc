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

#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>

namespace cplxfnc {

// ##################################################
// ##     Hurwitz Zeta function
// ##################################################

std::complex<double> zeta(std::complex<double> s, std::complex<double> a)
{
    return zeta(s, a, 1e-16, 5, false);
}

std::complex<double> zeta(std::complex<double> s, std::complex<double> a, double tol, unsigned int limit, bool verbose)
{
    std::complex<double> res;
    int status = zeta(s, a, &res, tol, limit, verbose);
    std::ostringstream oss;
    oss << "zeta did not converge for s=" << s << " and a=" << a << " , try dec. tol or inc. limit";
    return status ? throw std::runtime_error(oss.str()) : res;
}

int zeta(std::complex<double> s, std::complex<double> a, std::complex<double> * res, double tol,
         unsigned int limit, bool verbose)
{
    acb_t _z, _s, _a;
    acb_init(_z); acb_init(_s); acb_init(_a);
    acb_set_d_d(_s, s.real(), s.imag());
    acb_set_d_d(_a, a.real(), a.imag());

    unsigned int prec = 56;  /* this value was found to give error
                                below 1e-16 at least for s and a of the order of one  */
    unsigned int c = 1;

//    std::cout << "zeta: " << s << " , " << a << " ... ";

    double zeta_re, zeta_re_err, zeta_im, zeta_im_err, frac_re, frac_im;

    while (1) {
        acb_hurwitz_zeta(_z, _s, _a, prec);

        arb_ptr re = acb_realref(_z);
//        acb_print(_z);
//        std::cout << std::endl;
        zeta_re     = arf_get_d(arb_midref(re), ARF_RND_NEAR);
        zeta_re_err = mag_get_d(arb_radref(re));

        arb_ptr im = acb_imagref(_z);
        zeta_im     = arf_get_d(arb_midref(im), ARF_RND_NEAR);
        zeta_im_err = mag_get_d(arb_radref(im));

        if ((zeta_im == 0) && (zeta_im_err == 0)) {
            frac_im = 0;
        } else {
            frac_im = fabs(zeta_im_err) / fabs(zeta_im);
        }
        if ((zeta_re == 0) && (zeta_re_err == 0)) {
            frac_re = 0;
        } else {
            frac_re = fabs(zeta_re_err) / fabs(zeta_re);
        }

        if ( (frac_im < tol) && (frac_re < tol) ) {
            acb_clear(_z); acb_clear(_s); acb_clear(_a);
//            std::cout << std::complex<double>(zeta_re, zeta_im) << std::endl;
            *res = std::complex<double>(zeta_re, zeta_im);
            return 0;
        }
        prec *= 2;
        c += 1;
        if (c > limit) {
            if (verbose) {
                std::cerr << "zeta limit (" << limit << ") reached\n" <<
                std::setprecision(1) << std::fixed <<
                "zeta(s, a) with s=" << s << " and a=" << a << std::endl <<
                std::scientific << std::setprecision(2) <<
                "tol: " << tol << std::endl <<
                std::scientific << std::setprecision(16) <<
                "re : " << zeta_re << " +/- " << zeta_re_err << "  ->  frac: " << frac_re << std::endl <<
                "im : " << zeta_im << " +/- " << zeta_im_err << "  ->  frac: " << frac_im << std::endl;
            }
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
    if (status == 0) {
        return res;
    } else {
        std::ostringstream oss;
        if (status == -1) {
            oss << "gamma_inc did not converge for s=" << s << " and z=" << z << " , try dec. tol or inc. limit";
            throw std::runtime_error(oss.str());
        } else if (status == -2) {
            oss << "gamma_inc value error: if Re(s) < 0 then z must not be zero!";
            throw std::runtime_error(oss.str());
        } else {
            oss << "gamma_inc unknown error: error code: " << status;
            throw std::runtime_error(oss.str());
        }
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
    double res_re, res_re_err, res_im, res_im_err, frac_re, frac_im;

    while (1) {
        acb_hypgeom_gamma_upper(_res, _s, _z, 0, prec);   //

        arb_ptr re = acb_realref(_res);
        res_re     = arf_get_d(arb_midref(re), ARF_RND_NEAR);
        res_re_err = mag_get_d(arb_radref(re));

        arb_ptr im = acb_imagref(_res);
        res_im     = arf_get_d(arb_midref(im), ARF_RND_NEAR);
        res_im_err = mag_get_d(arb_radref(im));

        if ((res_im == 0) && (res_im_err == 0)) {
            frac_im = 0;
        } else {
            frac_im = fabs(res_im_err) / fabs(res_im);
        }
        if ((res_re == 0) && (res_re_err == 0)) {
            frac_re = 0;
        } else {
            frac_re = fabs(res_re_err) / fabs(res_re);
        }

        if ( (frac_im < tol) && (frac_re < tol) ) {
            acb_clear(_res); acb_clear(_s); acb_clear(_z);
            *res = std::complex<double>(res_re, res_im);
            return 0;
        }
        c += 1;
        if (c > limit) {
            if (verbose) {
                std::cerr << "ERROR: inc gamma limit (" << limit << ") reached\n" <<
                std::setprecision(1) << std::fixed <<
                "gamma(s, z) with s=" << s << " and z=" << z << std::endl <<
                "internal prec: " << prec << std::endl <<
                std::scientific << std::setprecision(2) <<
                "tol: " << tol << std::endl <<
                std::scientific << std::setprecision(16) <<
                "re : " << res_re << " +/- " << res_re_err << "  ->  frac: " << frac_re << std::endl <<
                "im : " << res_im << " +/- " << res_im_err << "  ->  frac: " << frac_im << std::endl;
            }
            return -1;
        } else {
            if (verbose) {
                std::cout << "gamma(s, z) with s=" << s << " and z=" << z << std::endl <<
                "internal prec: " << prec << std::endl <<
                std::scientific << std::setprecision(2) <<
                "tol: " << tol << std::endl <<
                std::scientific << std::setprecision(16) <<
                "re : " << res_re << " +/- " << res_re_err << "  ->  frac: " << frac_re << std::endl <<
                "im : " << res_im << " +/- " << res_im_err << "  ->  frac: " << frac_im << std::endl;
                std::cout << "new internal prec: " << prec << std::endl;

            }
        }
        prec *= 2;

    }
}

} /* namespace cplxfnc */