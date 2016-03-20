#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>

#include <stan/math.hpp>

// translation of Hornberg model to Stan, model from
// http://jjj.biochem.sun.ac.za/database/hornberg/index.html and
// described in paper DOI 10.1111/j.1432-1033.2004.04404.x
#include <test/unit/math/rev/arr/functor/hornberg.hpp>
#include <test/unit/math/rev/arr/functor/hornberg_jacobian.hpp>

// the rest of this test stays exactly the same
#include <test/unit/math/rev/arr/functor/integrate_ode_perf_test.cpp>

