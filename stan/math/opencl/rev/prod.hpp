#ifndef STAN_MATH_OPENCL_REV_PROD_HPP
#define STAN_MATH_OPENCL_REV_PROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/prod.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/opencl/rev/scalar_cl.hpp>

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
inline opencl::ScalarCl<var> prod(const var_value<T>& x) {
  return opencl::make_callback_scalar_cl(
      prod(value_of(x)), [x](const auto& res_adj, const auto& res_val) mutable {
        x.adj() += elt_divide(res_adj * res_val, x.val());
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
