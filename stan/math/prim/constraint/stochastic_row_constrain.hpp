#ifndef STAN_MATH_PRIM_CONSTRAINT_STOCHASTIC_ROW_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_STOCHASTIC_ROW_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/simplex_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a row stochastic matrix.
 *
 * The transform is defined using the inverse of the
 * isometric log ratio (ILR) transform
 *
 * @tparam Mat type of the Matrix
 * @param y Free Matrix input of dimensionality (N, K - 1).
 * @return Matrix with simplexes along the rows of dimensionality (N, K).
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_row_constrain(const Mat& y) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index N = y_ref.rows();
  plain_type_t<Mat> ret(N, y_ref.cols() + 1);
  for (Eigen::Index i = 0; i < N; ++i) {
    ret.row(i) = simplex_constrain(y_ref.row(i));
  }
  return ret;
}

/**
 * Return a row stochastic matrix.
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform
 *
 * @tparam Mat type of the matrix
 * @tparam Lp A scalar type for the lp argument. The scalar type of Mat should
 * be convertable to this.
 * @param y Free matrix input of dimensionality (N, K - 1).
 * @param lp Log probability reference to increment.
 * @return Matrix with simplexes along the rows of dimensionality (N, K).
 */
template <typename Mat, typename Lp,
          require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr,
          require_convertible_t<value_type_t<Mat>, Lp>* = nullptr>
inline plain_type_t<Mat> stochastic_row_constrain(const Mat& y, Lp& lp) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index N = y_ref.rows();
  plain_type_t<Mat> ret(N, y_ref.cols() + 1);
  for (Eigen::Index i = 0; i < N; ++i) {
    ret.row(i) = simplex_constrain(y_ref.row(i), lp);
  }
  return ret;
}

/**
 * Return a row stochastic matrix.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @param[in] y free vector with matrices of size (N, K - 1)
 * @return vector of matrices with simplex rows of dimensionality (N, K)
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_row_constrain(T&& y) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [](auto&& v) {
    return stochastic_row_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return a row stochastic matrix.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @tparam Lp Scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param[in] y free vector with matrices of size (N, K - 1)
 * @param[in, out] lp log density accumulator
 * @return vector of matrices with simplex rows of dimensionality (N, K)
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto stochastic_row_constrain(T&& y, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [&lp](auto&& v) {
    return stochastic_row_constrain(std::forward<decltype(v)>(v), lp);
  });
}

/**
 * Return a row stochastic matrix.
 * If the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Mat A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and dynamic columns, or a standard vector thereof
 * @tparam Lp A scalar type for the lp argument. The scalar type of Mat should
 * be convertable to this.
 * @param[in] y free matrix
 * @param[in, out] lp log density accumulator
 * @return Matrix with simplexes along the rows of dimensionality (N, K).
 */
template <bool Jacobian, typename Mat, typename Lp,
          require_convertible_t<return_type_t<Mat>, Lp>* = nullptr>
inline plain_type_t<Mat> stochastic_row_constrain(Mat&& y, Lp& lp) {
  if constexpr (Jacobian) {
    return stochastic_row_constrain(std::forward<Mat>(y), lp);
  } else {
    return stochastic_row_constrain(std::forward<Mat>(y));
  }
}

}  // namespace math
}  // namespace stan

#endif
