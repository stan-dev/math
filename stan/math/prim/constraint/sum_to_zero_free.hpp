#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
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

  y.coeffRef(N - 1) = -z_ref(N) * sqrt(N * (N + 1)) / N;

  value_type_t<Vec> sum_w(0);

  for (int i = N - 2; i >= 0; --i) {
    double n = static_cast<double>(i + 1);
    auto w = y(i + 1) / sqrt((n + 1) * (n + 2));
    sum_w += w;
    y.coeffRef(i) = (sum_w - z_ref(i + 1)) * sqrt(n * (n + 1)) / n;
  }

  return y;
}

/**
 * Overload of `sum_to_zero_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param z The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto sum_to_zero_free(const T& z) {
  return apply_vector_unary<T>::apply(
      z, [](auto&& v) { return sum_to_zero_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
