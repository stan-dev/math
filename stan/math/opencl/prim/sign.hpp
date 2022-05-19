#ifndef STAN_MATH_OPENCL_PRIM_SIGN_HPP
#define STAN_MATH_OPENCL_PRIM_SIGN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Returns signs of the arguments.
 * @tparam Matcl type of the argument (`matrix_cl` or kernel generator
 * expression)
 * @param x the argument
 * @return sign of `x`
 */
template <typename Matcl,
          require_all_kernel_expressions_and_none_scalar_t<Matcl>* = nullptr>
auto sign(const Matcl& x) {
  return select(x == 0, 0, select(x < 0, -1, 1));
}

}  // namespace math
}  // namespace stan
#endif
#endif
