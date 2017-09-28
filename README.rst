cplxfnc -- special functions with complex arguments based on arb library
========================================================================


This library allows to evaluate special function with complex arguments up to double precision
The arb library (http://arblib.org) is used behind the scenes which allows for correct error estimation. The internal precision
is increased until the result has converged to double precision.

wrapped functions:
  * Hurwitz zeta function