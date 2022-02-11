#ifndef STAN_MATH_OPENCL_PRIM_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_OPENCL_PRIM_COLUMNS_DOT_SELF_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of the matrix
 *
 * @param a Matrix.
 * @return Row vector containing the dot product of each column of the matrix
 * with itself.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto columns_dot_self(const T& a) {
  if (size_zero(a)) {
    return plain_type_t<T>(constant(0.0, 1, a.cols()));
  }
  plain_type_t<T> res = colwise_sum(square(a));
  while (res.rows() > 1) {
    res = colwise_sum(res).eval();
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
