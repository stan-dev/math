#ifndef STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_COLUMN_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_COLUMN_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/simplex_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a column stochastic matrix.
 *
 * The transform is defined using the inverse of the
 * isometric log ratio (ILR) transform
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
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform
 *
 * @tparam Mat type of the Matrix
 * @tparam Lp A scalar type for the lp argument. The scalar type of Mat should
 * be convertable to this.
 * @param y Free Matrix input of dimensionality (K - 1, M)
 * @param lp Log probability reference to increment.
 * @return Matrix with stochastic columns of dimensionality (K, M)
 */
template <typename Mat, typename Lp,
          require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr,
          require_convertible_t<value_type_t<Mat>, Lp>* = nullptr>
inline plain_type_t<Mat> stochastic_column_constrain(const Mat& y, Lp& lp) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index M = y_ref.cols();
  plain_type_t<Mat> ret(y_ref.rows() + 1, M);
  for (Eigen::Index i = 0; i < M; ++i) {
    ret.col(i) = simplex_constrain(y_ref.col(i), lp);
  }
  return ret;
}
/**
 * Return a vector of column stochastic matrices.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @param[in] y free vector
 * @return Standard vector containing matrices with simplex columns of
 * dimensionality (K, M).
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_column_constrain(T&& y) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [](auto&& v) {
    return stochastic_column_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return a vector of column stochastic matrices.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @tparam Lp Scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return Standard vector containing matrices with simplex columns of
 * dimensionality (K, M).
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto stochastic_column_constrain(T&& y, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [&lp](auto&& v) {
    return stochastic_column_constrain(std::forward<decltype(v)>(v), lp);
  });
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
 * @tparam Mat A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and dynamic columns, or a standard vector thereof
 * @tparam Lp A scalar type for the lp argument. The scalar type of Mat should
 * be convertable to this.
 * @param y Free Matrix input of dimensionality (K - 1, M).
 * @param[in, out] lp log density accumulator
 * @return Matrix with simplex columns of dimensionality (K, M).
 */
template <bool Jacobian, typename Mat, typename Lp,
          require_convertible_t<return_type_t<Mat>, Lp>* = nullptr>
inline plain_type_t<Mat> stochastic_column_constrain(Mat&& y, Lp& lp) {
  if constexpr (Jacobian) {
    return stochastic_column_constrain(std::forward<Mat>(y), lp);
  } else {
    return stochastic_column_constrain(std::forward<Mat>(y));
  }
}

}  // namespace math
}  // namespace stan

#endif
