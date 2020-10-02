#ifndef STAN_MATH_OPENCL_REV_SUM_HPP
#define STAN_MATH_OPENCL_REV_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Returns the sum of the coefficients of the specified
 * matrix on the OpenCL device.
 *
 * @param x Specified var_value containing a matrix.
 * @return Sum of coefficients of matrix.
 */
inline var sum(const var_value<matrix_cl<double>>& x) {
  var res = sum(value_of(x));
  reverse_pass_callback([res, x]() mutable { x.adj() = x.adj() + res.adj(); });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
