#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of the
 * the inverse of the isometric log ratio transform (ILR).
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
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Zero-sum vector of dimensionality K.
 */
template <typename Vec, require_eigen_col_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y) {
  const auto N = y.size();

  plain_type_t<Vec> z = Eigen::VectorXd::Zero(N + 1);
  if (unlikely(N == 0)) {
    return z;
  }

  auto&& y_ref = to_ref(y);

  value_type_t<Vec> sum_w(0);
  for (int i = N; i > 0; --i) {
    double n = static_cast<double>(i);
    auto w = y_ref.coeff(i - 1) * inv_sqrt(n * (n + 1));
    sum_w += w;

    z.coeffRef(i - 1) += sum_w;
    z.coeffRef(i) -= w * n;
  }

  return z;
}

/**
 * Return a matrix that sums to zero over both the rows
 * and columns corresponding to the free matrix x.
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Mat type of the matrix
 * @param x Free matrix input of dimensionality (N - 1, M - 1).
 * @return Zero-sum matrix of dimensionality (N, M).
 */
template <typename Mat, require_eigen_matrix_dynamic_t<Mat>* = nullptr,
          require_not_st_var<Mat>* = nullptr>
inline plain_type_t<Mat> sum_to_zero_constrain(const Mat& x) {
  const auto N = x.rows();
  const auto M = x.cols();

  plain_type_t<Mat> Z = Eigen::MatrixXd::Zero(N + 1, M + 1);
  if (unlikely(N == 0 || M == 0)) {
    return Z;
  }
  auto&& x_ref = to_ref(x);

  Eigen::Matrix<value_type_t<Mat>, -1, 1> beta = Eigen::VectorXd::Zero(N);

  for (int j = M - 1; j >= 0; --j) {
    value_type_t<Mat> ax_previous(0);

    double a_j = inv_sqrt((j + 1.0) * (j + 2.0));
    double b_j = (j + 1.0) * a_j;

    for (int i = N - 1; i >= 0; --i) {
      double a_i = inv_sqrt((i + 1.0) * (i + 2.0));
      double b_i = (i + 1.0) * a_i;

      auto b_i_x = b_i * x_ref.coeff(i, j) - ax_previous;

      Z.coeffRef(i, j) = (b_j * b_i_x) - beta.coeff(i);
      beta.coeffRef(i) += a_j * b_i_x;

      Z.coeffRef(N, j) -= Z.coeff(i, j);
      Z.coeffRef(i, M) -= Z.coeff(i, j);

      ax_previous += a_i * x_ref.coeff(i, j);
    }
    Z.coeffRef(N, M) -= Z.coeff(N, j);
  }

  return Z;
}

/**
 * Return a vector or matrix with sum zero corresponding to the specified
 * free input.
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam T type of the input, either a vector or a matrix
 * @tparam Lp unused
 * @param y Free vector or matrix
 * @param lp unused
 * @return Zero-sum vector or matrix which is one larger in each dimension
 */
template <typename T, typename Lp, require_not_st_var<T>* = nullptr>
inline plain_type_t<T> sum_to_zero_constrain(T&& y, Lp& lp) {
  return sum_to_zero_constrain(std::forward<T>(y));
}

/**
 * Return a vector or matrix with sum zero corresponding to the specified
 * free input.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type that is either a vector or a
 * matrix
 * @param[in] y free vector or matrix
 * @return Zero-sum vectors or matrices which are one larger in each dimension
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(T&& y) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [](auto&& v) {
    return sum_to_zero_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type that is either a vector or a
 * matrix
 * @tparam Lp unused
 * @param[in] y free vector or matrix
 * @param[in, out] lp unused
 * @return Zero-sum vectors or matrices which are one larger in each dimension
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto sum_to_zero_constrain(T&& y, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [](auto&& v) {
    return sum_to_zero_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return a vector or matrix with sum zero corresponding to the specified
 * free input.
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Jacobian unused
 * @tparam T type of the input
 * @tparam Lp unused
 * @param[in] y free vector or matrix
 * @param[in, out] lp unused
 * @return Zero-sum vector or matrix which is one larger in each dimension
 */
template <bool Jacobian, typename T, typename Lp>
inline plain_type_t<T> sum_to_zero_constrain(T&& y, Lp& lp) {
  return sum_to_zero_constrain(std::forward<T>(y));
}

}  // namespace math
}  // namespace stan

#endif
