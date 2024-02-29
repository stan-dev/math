#ifndef STAN_MATH_PRIM_FUN_STOCHASTIC_COLUMN_FREE_HPP
#define STAN_MATH_PRIM_FUN_STOCHASTIC_COLUMN_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/simplex_free.hpp>

namespace stan {
namespace math {

/**
 * Return an unconstrained matrix that when transformed produces
 * the specified columnwise stochastic matrix.  It applies to a stochastic
 * matrix of dimensionality (N, K) and produces an unconstrained vector of
 * dimensionality (N - 1, K).
 *
 * @tparam Mat type of the Matrix
 * @param y Columnwise stochastic matrix input of dimensionality (N, K)
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_column_free(const Mat& y) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index M = y_ref.cols();
  plain_type_t<Mat> ret(y_ref.rows() - 1, M);
  for (Eigen::Index i = 0; i < M; ++i) {
    ret.col(i) = simplex_free(y_ref.col(i));
  }
  return ret;
}

/**
 * Overload that untransforms each Eigen matrix in a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic rows
 * @param[in] y vector of columnwise stochastic matrix of size (N, K)
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_column_free(const T& y) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return stochastic_column_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
