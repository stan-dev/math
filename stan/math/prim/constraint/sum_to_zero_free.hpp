#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an unconstrained vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements.
 *
 * @tparam ColVec a column vector type
 * @param x Vector of length K.
 * @return Free vector of length (K-1).
 * @throw std::domain_error if x does not sum to zero
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_free(const Vec& x) {
  const auto& x_ref = to_ref(x);
  check_sum_to_zero("stan::math::sum_to_zero_free", "sum_to_zero variable",
                    x_ref);
  if (x_ref.size() == 0) {
    return plain_type_t<Vec>(0);
  }
  return x_ref.head(x_ref.size() - 1);
}

/**
 * Overload of `sum_to_zero_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto sum_to_zero_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return sum_to_zero_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
