#ifndef STAN_MATH_OPENCL_PRIM_SQUARED_DISTANCE_HPP
#define STAN_MATH_OPENCL_PRIM_SQUARED_DISTANCE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 *
 * @tparam T_a type of the first kernel generator expression
 * @tparam T_b type of the second kernel generator expression
 *
 * @param a first kernel generator expression
 * @param b second kernel generator expression
 * @return squared distance.
 * @throw std::domain_error If the matrices are not the same
 * size
 */
template <typename T_a, typename T_b,
          require_all_kernel_expressions_and_none_scalar_t<T_a, T_b>* = nullptr>
inline auto squared_distance(const T_a& a, const T_b& b) {
  const char* function = "squared_distance (OpenCL)";
  check_vector(function, "a", a);
  check_vector(function, "b", b);
  return sum(square(a - b));
}

}  // namespace math
}  // namespace stan
#endif
#endif
