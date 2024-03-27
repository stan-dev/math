#ifndef STAN_MATH_PRIM_CONSTRAINT_STOCHASTIC_ROW_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_STOCHASTIC_ROW_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/simplex_free.hpp>

namespace stan {
namespace math {

/**
 * Return an unconstrained matrix that when transformed produces
 * the specified simplex matrix.  It applies to a simplex of dimensionality
 * (N, K) and produces an unconstrained vector of dimensionality (N, K - 1).
 *
 * @tparam Mat type of the Matrix
 * @param y Rowwise simplex Matrix input of dimensionality (N, K)
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> stochastic_row_free(const Mat& y) {
  auto&& y_ref = to_ref(y);
  const Eigen::Index N = y_ref.rows();
  plain_type_t<Mat> ret(N, y_ref.cols() - 1);
  for (Eigen::Index i = 0; i < N; ++i) {
    ret.row(i) = simplex_free(y_ref.row(i));
  }
  return ret;
}

/**
 * Overload that untransforms each Eigen matrix in a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic rows
 * @param[in] y vector of rowwise simplex matrices each of size (N, K)
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_row_free(const T& y) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return stochastic_row_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
