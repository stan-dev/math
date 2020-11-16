#ifndef STAN_MATH_OPENCL_PRIM_SIGN_HPP
#define STAN_MATH_OPENCL_PRIM_SIGN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Returns signs of the arguments.
 * @tparam T type of the argument (`matrix_cl` or kernel generator expression)
 * @param x the argument
 * @return sign of `x`
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
auto sign(const T& x) {
  return select(x == 0, 0, select(x < 0, -1, 1));
}

}  // namespace math
}  // namespace stan
#endif
#endif
