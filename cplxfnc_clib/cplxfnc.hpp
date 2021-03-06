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

#ifndef CPLXFNC_H
#define CPLXFNC_H

#include <complex>

#define ZETA_DEFAULT_INIT_PREC 56
#define GAMMA_INC_DEFAULT_INIT_PREC 75
#define U_ASYMP_DEFAULT_INIT_PREC 56

extern "C" {
  void libcplxfnc_is_present(void);
}

namespace cplxfnc {

std::complex<double> zeta(std::complex<double> s, std::complex<double> a);
std::complex<double> zeta(std::complex<double> s, std::complex<double> a, double tol, 
        unsigned int limit, bool verbose,
        unsigned int init_prec=ZETA_DEFAULT_INIT_PREC);
int zeta(std::complex<double> s, std::complex<double> a, std::complex<double> *res, double tol,
        unsigned int limit, bool verbose, 
        unsigned int init_prec=ZETA_DEFAULT_INIT_PREC);

std::complex<double> gamma_inc(std::complex<double> s, std::complex<double> z);
std::complex<double> gamma_inc(std::complex<double> s, std::complex<double> z, double tol,
        unsigned int limit, bool verbose,
        unsigned int init_prec=GAMMA_INC_DEFAULT_INIT_PREC);
int gamma_inc(std::complex<double> s, std::complex<double> z, std::complex<double> * res, double tol,
        unsigned int limit, bool verbose, 
        unsigned int init_prec=GAMMA_INC_DEFAULT_INIT_PREC);

std::complex<double> u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z);              
std::complex<double> u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z, double tol,
        unsigned int limit, bool verbose, 
        unsigned int init_prec=U_ASYMP_DEFAULT_INIT_PREC);
int u_asymp(std::complex<double> a, std::complex<double> b, std::complex<double> z, std::complex<double> * res, double tol,
        unsigned int limit, bool verbose, 
        unsigned int init_prec=U_ASYMP_DEFAULT_INIT_PREC);

}

#endif
