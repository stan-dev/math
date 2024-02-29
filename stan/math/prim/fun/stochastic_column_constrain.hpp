#ifndef STAN_MATH_PRIM_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/simplex_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a column stochastic matrix.
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam Mat type of the Matrix
 * @param y Free Matrix input of dimensionality (K - 1, M)
 * @return Matrix with simplex columns of dimensionality (K, M)
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_column_constrain(const Mat& y) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index M = y_ref.cols();
  plain_type_t<Mat> ret(y_ref.rows() + 1, M);
  for (Eigen::Index i = 0; i < M; ++i) {
    ret.col(i) = simplex_constrain(y_ref.col(i));
  }
  return ret;
}

/**
 * Return a column stochastic matrix
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam Mat type of the Matrix
 * @param y Free Matrix input of dimensionality (K - 1, M)
 * @param lp Log probability reference to increment.
 * @return Matrix with stochastic columns of dimensionality (K, M)
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_column_constrain(const Mat& y,
                                                     value_type_t<Mat>& lp) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index M = y_ref.cols();
  plain_type_t<Mat> ret(y_ref.rows() + 1, M);
  for (Eigen::Index i = 0; i < M; ++i) {
    ret.col(i) = simplex_constrain(y_ref.col(i), lp);
  }
  return ret;
}

/**
 * Return a column stochastic matrix. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Mat type of the Matrix
 * @param y Free Matrix input of dimensionality (K - 1, M).
 * @param[in, out] lp log density accumulator
 * @return Matrix with simplex columns of dimensionality (K, M).
 */
template <bool Jacobian, typename Mat, require_not_std_vector_t<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_column_constrain(const Mat& y,
                                                     return_type_t<Mat>& lp) {
  if (Jacobian) {
    return stochastic_column_constrain(y, lp);
  } else {
    return stochastic_column_constrain(y);
  }
}

/**
 * Return a vector of column stochastic matrices. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return Standard vector containing matrices with simplex columns of
 * dimensionality (K, M).
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_column_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(y, [&lp](auto&& v) {
    return stochastic_column_constrain<Jacobian>(v, lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
