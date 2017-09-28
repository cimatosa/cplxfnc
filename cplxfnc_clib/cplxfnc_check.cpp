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

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

int check_values(){
    std::cout << "check values ... ";

    std::complex<double> res, s, a, res_check;
    double d;
    double tol = 1e-16;
    const std::complex<double> I(0, 1);

    unsigned int limit = 1;

    try {
        res = cplxfnc::zeta(2., 1. - I*1234., tol, limit, true);
    } catch  (const std::runtime_error& e) {
        std::cout << "\nI guess you have arb version prior to 2.9. installed" <<
                     "\nupdating yields a slight efficiency improvement!\n";
        limit = 2;
    }

    double data [18][6] = { {2. , 0., 1.,    0. , M_PI*M_PI/6, 0},
                            {2. , 0., 1.,    1. ,  4.63000096622763786e-1, -7.94233542759318865e-1},
                            {2. , 0., 0.,    1. , -5.36999903377236213e-1, -7.94233542759318865e-1},
                            {2. , 0., 1.,   10. ,  4.99999999999999999e-3, -9.98329975849301549e-2},
                            {2. , 0., 1.,-1234. ,  3.28352014373937781e-7,  8.10372682779022825e-4},
                            {1.2, 0., 1.,    0. ,  5.59158244117775188   ,  0},
                            {1.2, 0., 1.,    1. ,  4.79688290618251932   , -1.04259444248270157},
                            {1.2, 0., 0.,    1. ,  4.48786591180757196   , -1.99365095877785516},
                            {1.2, 0., 1.,   10. ,  3.00952851308319147   , -9.44683699633929863e-1},
                            {1.2, 0., 1.,-1234. ,  1.14531441635402694   ,  3.72032603245348282e-1},
                            {1.2, 1., 1.,    0. ,  7.89008276053315959e-1, -8.904382488200783911e-1},
                            {1.2, 1., 1.,    1. , -3.74965169663485969e-1, -2.732104072569482620},
                            {1.2, 1., 0.,    1. , -1.86148443143806753   , -7.307139932227024369},
                            {1.2, 1., 1.,   10. , -1.89011737420273344   ,  2.105544806767073391},
                            {1.2, 1., 1.,-1234. , -1.56061472179186142e-2, -4.656891594125561926e-2},
                            {1.2, 0., 1.,  1.e6 ,  3.00038056734908233e-1, -9.748824108122723291e-2},
                            {1.2, 0., 1., 1.e12 ,  1.89311209369371249e-2, -6.151094064175929217e-3},
                            {1.2, 0., 1., 1.e24 ,  7.53661499160987508e-5, -2.448794653698247891e-5} };

    for (int i = 0; i < 18; i++){
        s         = data[i][0] + I*data[i][1];
        a         = data[i][2] + I*data[i][3];
        res_check = data[i][4] + I*data[i][5];

        res = cplxfnc::zeta(s, a, tol, limit, true);
        d = fabs(res.real() - res_check.real());
        if (d > tol){
            std::cout << "\nERROR (real part)\n" <<
                 std::scientific << std::setprecision(2) <<
                "diff real part: " << d << " < tol (" << tol << ")\n" <<
                 std::setprecision(1) << std::fixed <<
                "zeta(s, a) with s=" << s << " and a=" << a << std::endl << std::scientific << std::setprecision(16) <<
                "returned      : " << res.real() << std::endl <<
                "but should be : " << res_check.real() << std::endl;
            return -1;
        }
        d = fabs(res.imag() - res_check.imag());
        if (d > tol){
            std::cout << "\nERROR (imag part)\n" <<
                 std::scientific << std::setprecision(2) <<
                "diff imag part: " << d << " < tol (" << tol << ")\n" <<
                 std::setprecision(1) << std::fixed <<
                "zeta(s, a) with s=" << s << " and a=" << a << std::endl << std::scientific << std::setprecision(16) <<
                "returned      : " << res.imag() << std::endl <<
                "but should be : " << res_check.imag() << std::endl;
            return -1;
        }

    }
    std::cout << "done\n";
    return 0;
}

int check_call_error()
{
    std::cout << "check call err ... ";

    std::complex<double> res, s, a, res_check;
    double tol = 1e-16;
    const std::complex<double> I(0, 1);

    s = 1.00000004 + I*10.;
    a = 1.e6 + I*1.e6;

    int status = cplxfnc::zeta(s, a, &res, tol, 1, false);
    if (status == 0) {
        std::cout << "\nERROR (limit should have been reached)\n" <<
        "expect zeta to return other than 0 return code!" << std::endl;
        return -1;
    }


    bool got_exception = false;
    try {
        res = cplxfnc::zeta(s, a, tol, 1, false);
    } catch (const std::runtime_error& e) {
        got_exception = true;
    }

    if (got_exception == false) {
        std::cout << "\nERROR (limit should have been reached)\n" <<
        "expect zeta to raise runtime_error exception!" << std::endl;
        return -1;
    }

    std::cout << "done\n";
    return 0;
}

int check_call_overloads()
{
    std::cout << "check call overloads ... ";
    std::complex<double> res;
    cplxfnc::zeta(2, 1, &res, 1e-16, 1, false);
    cplxfnc::zeta(2, 1, 1e-16, 1, false);
    cplxfnc::zeta(2, 1);
    std::cout << "done\n";
    return 0;
}

int main(){
    std::cout << "\nrun tests for cplxfnc library\n";

    if (check_values()) return -1;
    if (check_call_error()) return -1;
    if (check_call_overloads()) return -1;

    return 0;
}
