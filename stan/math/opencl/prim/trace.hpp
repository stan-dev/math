#ifndef STAN_MATH_OPENCL_PRIM_TRACE_HPP
#define STAN_MATH_OPENCL_PRIM_TRACE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>

namespace stan {
namespace math {

/**
 * Calculates trace (sum of diagonal) of given kernel generator expression.
 * @tparam T type of the expression
 * @param m expression to calculate trace of
 * @return trace of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
value_type_t<T> trace(const T& m) {
  return sum(diagonal(m));
}

}  // namespace math
}  // namespace stan

#endif
#endif
