#ifndef STAN_MATH_REV_SUNDIALS_CHECK_HPP
#define STAN_MATH_REV_SUNDIALS_CHECK_HPP

#include <stan/math/prim/err/check_flag_sundials.hpp>

#define CHECK_SUNDIALS_CALL(call) check_flag_sundials(call, #call)

#endif
