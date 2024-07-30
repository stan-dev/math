#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Zero-sum vector of dimensionality K.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y) {
  using T = value_type_t<Vec>;

  int Km1 = y.size();
  plain_type_t<Vec> x(Km1 + 1);
  // copy the first Km1 elements
  x.head(Km1) = y;
  // set the last element to -sum(y)
  x.coeffRef(Km1) = -sum(y);
  return x;
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements. This is a linear transform, with no
 * Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @param lp unused
 * @return Zero-sum vector of dimensionality K.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y,
                                               value_type_t<Vec>& lp) {
  return sum_to_zero_constrain(y);
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements. This is a linear transform, with no
 * Jacobian.
 *
 * @tparam Jacobian unused
 * @tparam Vec A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param[in] y free vector
 * @param[in, out] lp unused
 * @return Zero-sum vector of dimensionality one greater than `y`
 */
template <bool Jacobian, typename Vec, require_not_std_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y,
                                               return_type_t<Vec>& lp) {
  return sum_to_zero_constrain(y);
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements. This is a linear transform, with no
 * Jacobian.
 *
 * @tparam Jacobian unused
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param[in] y free vector
 * @param[in, out] lp unused
 * @return Zero-sum vectors of dimensionality one greater than `y`
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return sum_to_zero_constrain(v); });
}

}  // namespace math
}  // namespace stan

#endif
