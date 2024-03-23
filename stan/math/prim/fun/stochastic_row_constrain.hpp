#ifndef STAN_MATH_PRIM_FUN_STOCHASTIC_ROW_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_STOCHASTIC_ROW_CONSTRAIN_HPP

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
 * Return a row stochastic matrix.
 *
 * The transform is based on a centered stick-breaking process.
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
  int Km1 = y_ref.cols();
  plain_type_t<Mat> x(N, Km1 + 1);
  using eigen_arr = Eigen::Array<scalar_type_t<Mat>, -1, 1>;
  eigen_arr stick_len = eigen_arr::Constant(N, 1.0);
  for (Eigen::Index k = 0; k < Km1; ++k) {
    auto z_k = inv_logit(y_ref.array().col(k) - log(Km1 - k));
    x.array().col(k) = stick_len * z_k;
    stick_len -= x.array().col(k);
  }
  x.array().col(Km1) = stick_len;
  return x;
}

/**
 * Return a row stochastic matrix.
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam Mat type of the matrix
 * @param y Free matrix input of dimensionality (N, K - 1).
 * @param lp Log probability reference to increment.
 * @return Matrix with simplexes along the rows of dimensionality (N, K).
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_row_constrain(const Mat& y,
                                                  value_type_t<Mat>& lp) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index N = y_ref.rows();
  Eigen::Index Km1 = y_ref.cols();
  plain_type_t<Mat> x(N, Km1 + 1);
  Eigen::Array<scalar_type_t<Mat>, -1, 1> stick_len
      = Eigen::Array<scalar_type_t<Mat>, -1, 1>::Constant(N, 1.0);
  for (Eigen::Index k = 0; k < Km1; ++k) {
    const auto eq_share = -log(Km1 - k);  // = logit(1.0/(Km1 + 1 - k));
    auto adj_y_k = (y_ref.array().col(k) + eq_share).eval();
    auto z_k = inv_logit(adj_y_k);
    x.array().col(k) = stick_len * z_k;
    lp += -sum(log1p_exp(adj_y_k)) - sum(log1p_exp(-adj_y_k))
          + sum(log(stick_len));
    stick_len -= x.array().col(k);  // equivalently *= (1 - z_k);
  }
  x.col(Km1).array() = stick_len;
  return x;
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
 *  and dynamic columns
 * @param[in] y free matrix
 * @param[in, out] lp log density accumulator
 * @return Matrix with simplexes along the rows of dimensionality (N, K).
 */
template <bool Jacobian, typename Mat, require_not_std_vector_t<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_row_constrain(const Mat& y,
                                                  return_type_t<Mat>& lp) {
  if (Jacobian) {
    return stochastic_row_constrain(y, lp);
  } else {
    return stochastic_row_constrain(y);
  }
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
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @param[in] y free vector with matrices of size (N, K - 1)
 * @param[in, out] lp log density accumulator
 * @return vector of matrices with simplex rows of dimensionality (N, K)
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_row_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      y, [&lp](auto&& v) { return stochastic_row_constrain<Jacobian>(v, lp); });
}

}  // namespace math
}  // namespace stan

#endif
