cplxfnc -- special functions with complex arguments based on arb library
========================================================================


This library allows to evaluate special function with complex arguments up to double precision.
The arb library (http://arblib.org) is used behind the scenes which allows for correct error estimation. The internal precision
is increased until the result has converged to double precision.

This package provides an interface for c++ as well as python.

## wrapped functions:

* Hurwitz zeta function

more to follow
  
## building the c++ code:  

The arb library needs to be installed. See http://arblib.org/setup.html for instructions. 
On Debian you may install `libflint-arb-dev` (included since Debian 9 or Ubuntu 17.04).
    
Then simply run the make script by typing
       
    cd ./cplxfnc_clib
    ./configure
    make
    
To install type `make install` (root permission). You also may run sanity checks: `make check`.