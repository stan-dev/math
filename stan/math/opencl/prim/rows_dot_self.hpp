#ifndef STAN_MATH_OPENCL_PRIM_ROWS_DOT_SELF_HPP
#define STAN_MATH_OPENCL_PRIM_ROWS_DOT_SELF_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 *
 * @tparam T_a type of the matrix
 *
 * @param a matrix
 * @return Dot product of the vectors.
 */
template <typename T_a,
          require_all_kernel_expressions_and_none_scalar_t<T_a>* = nullptr>
inline auto rows_dot_self(T_a&& a) {
  return rowwise_sum(square(std::forward<T_a>(a)));
}

}  // namespace math
}  // namespace stan

#endif
#endif
