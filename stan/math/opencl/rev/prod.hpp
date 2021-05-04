#ifndef STAN_MATH_OPENCL_REV_PROD_HPP
#define STAN_MATH_OPENCL_REV_PROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/prod.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Returns the prod of the coefficients of the specified
 * matrix on the OpenCL device.
 *
 * @param x Specified var_value containing a matrix.
 * @return prod of coefficients of matrix.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var prod(const var_value<T>& x) {
  return make_callback_var(prod(value_of(x)), [x](vari& res) mutable {
    x.adj() += elt_divide(res.adj() * res.val(), x.val());
  });
}

}  // namespace math
}  // namespace stan

#endif
#endif
