#ifndef STAN_MATH_OPENCL_PRIM_DOT_PRODUCT_HPP
#define STAN_MATH_OPENCL_PRIM_DOT_PRODUCT_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of the specified vectors.
 *
 * @tparam T_a type of the first vector
 * @tparam T_b type of the second vector
 * @param a First vector.
 * @param b Second vector.
 * @return Dot product of the vectors.
 * @throw std::invalid_argument If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T_a, typename T_b,
          require_all_kernel_expressions_and_none_scalar_t<T_a, T_b>* = nullptr>
inline auto dot_product(const T_a& a, const T_b& b) {
  const char* function = "dot_product(OpenCL)";
  check_vector(function, "a", a);
  check_vector(function, "b", b);
  return sum(elt_multiply(a, b));
}

}  // namespace math
}  // namespace stan

#endif
#endif
