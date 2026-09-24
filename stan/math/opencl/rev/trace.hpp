#ifndef STAN_MATH_OPENCL_REV_TRACE_HPP
#define STAN_MATH_OPENCL_REV_TRACE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/prim/trace.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/opencl/rev/scalar_cl.hpp>

namespace stan {
namespace math {

/**
 * Returns the trace (sum of diagonal) of the specified
 * matrix on the OpenCL device.
 *
 * @param x Specified var_value containing a matrix.
 * @return trace of matrix
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline opencl::ScalarCl<var> trace(const var_value<T>& x) {
  return opencl::make_callback_scalar_cl(
      trace(value_of(x)), [x](const auto& res_adj, const auto&) mutable {
        diagonal(x.adj()) += res_adj;
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
