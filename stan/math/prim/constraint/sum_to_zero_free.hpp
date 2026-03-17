#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an unconstrained vector.
 *
 * The sum-to-zero transform is defined using a modified version of the
 * isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * @tparam ColVec a column vector type
 * @param z Vector of length K.
 * @return Free vector of length (K-1).
 * @throw std::domain_error if z does not sum to zero
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_free(const Vec& z) {
  const auto& z_ref = to_ref(z);
  check_sum_to_zero("stan::math::sum_to_zero_free", "sum_to_zero variable",
                    z_ref);

  const auto N = z.size() - 1;

  plain_type_t<Vec> y = Eigen::VectorXd::Zero(N);
  if (unlikely(N == 0)) {
    return y;
  }

  y.coeffRef(N - 1) = -z_ref.coeff(N) * sqrt(N * (N + 1)) / N;

  value_type_t<Vec> sum_w(0);

  for (int i = N - 2; i >= 0; --i) {
    double n = static_cast<double>(i + 1);
    auto w = y.coeff(i + 1) / sqrt((n + 1) * (n + 2));
    sum_w += w;
    y.coeffRef(i) = (sum_w - z_ref.coeff(i + 1)) * sqrt(n * (n + 1)) / n;
  }

  return y;
}

/**
 * Return an unconstrained matrix.
 *
 * @tparam Mat a column vector type
 * @param z Matrix of size (N, M)
 * @return Free matrix of length (N - 1, M - 1)
 * @throw std::domain_error if z does not sum to zero
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr>
inline plain_type_t<Mat> sum_to_zero_free(const Mat& z) {
  const auto& z_ref = to_ref(z);
  check_sum_to_zero("stan::math::sum_to_zero_free", "sum_to_zero variable",
                    z_ref);

  const auto N = z_ref.rows() - 1;
  const auto M = z_ref.cols() - 1;

  plain_type_t<Mat> x = Eigen::MatrixXd::Zero(N, M);
  if (unlikely(N == 0 || M == 0)) {
    return x;
  }

  Eigen::Matrix<value_type_t<Mat>, -1, 1> beta = Eigen::VectorXd::Zero(N);

  for (int j = M - 1; j >= 0; --j) {
    value_type_t<Mat> ax_previous(0);

    double a_j = inv_sqrt((j + 1.0) * (j + 2.0));
    double b_j = (j + 1.0) * a_j;

    for (int i = N - 1; i >= 0; --i) {
      double a_i = inv_sqrt((i + 1.0) * (i + 2.0));
      double b_i = (i + 1.0) * a_i;

      auto alpha_plus_beta = z_ref.coeff(i, j) + beta.coeff(i);

      x.coeffRef(i, j) = (alpha_plus_beta + b_j * ax_previous) / (b_j * b_i);
      beta.coeffRef(i) += a_j * (b_i * x.coeff(i, j) - ax_previous);
      ax_previous += a_i * x.coeff(i, j);
    }
  }

  return x;
}

/**
 * Overload of `sum_to_zero_free()` to untransform each Eigen type
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`
 * @param z The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto sum_to_zero_free(T&& z) {
  return apply_vector_unary<T>::apply(std::forward<T>(z), [](auto&& v) {
    return sum_to_zero_free(std::forward<decltype(v)>(v));
  });
}

}  // namespace math
}  // namespace stan

#endif
