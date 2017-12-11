cplxfnc -- special functions with complex arguments based on arb library
========================================================================

![Travis CI][https://travis-ci.org/cimatosa/cplxfnc.svg?branch=master] ![codecov][https://codecov.io/gh/cimatosa/cplxfnc/branch/master/graph/badge.svg]

This library allows to evaluate special function with complex arguments up to double precision.
The arb library (http://arblib.org) is used behind the scenes which allows for correct error estimation. The internal precision
is increased until the result has converged to double precision.

This package provides an interface for c++ as well as python.

## wrapped functions:

* Hurwitz zeta function (zeta)
* upper incomplete gamma function (gamma_inc)

more to follow

## building the c++ code:  

The arb library needs to be installed. See http://arblib.org/setup.html for instructions. 
On Debian you may install `libflint-arb-dev` (included since Debian 9 or Ubuntu 17.04).
    
Then simply run the make script by typing
       
    cd ./cplxfnc_clib
    ./configure
    make
    
To install type `make install` (root permission). You also may run sanity checks: `make check`.

## install python extention

Running the `setup.py` script as follows builds and installs the the python package `cplxfnc`.
Note, the build process needs the shared libraries for flint (`libflint`) and arb (`libarb`, or when installed from
linux package `libflint-arb`).

    python setup.py build
    python setup.py install

To check if things are working correctly, make sure `pytest` and `mpmath` are installed and type

    py.test
    
## example

Once installed, the usage it straight forward.

    >>> import cplxfnc
    
    >>> cplxfnc.zeta(s=1+1j, a=3-5j)
    (-0.3269595185571998+0.04885844807914104j)
