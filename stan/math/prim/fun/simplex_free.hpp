#ifndef STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP
#define STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an unconstrained vector that when transformed produces
 * the specified simplex.  It applies to a simplex of dimensionality
 * K and produces an unconstrained vector of dimensionality (K-1).
 *
 * <p>The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam ColVec type of the simplex (must be a column vector)
 * @param x Simplex of dimensionality K.
 * @return Free vector of dimensionality (K-1) that transforms to
 * the simplex.
 * @throw std::domain_error if x is not a valid simplex
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> simplex_free(const Vec& x) {
  using std::log;
  using T = value_type_t<Vec>;

  const auto& x_ref = to_ref(x);
  check_simplex("stan::math::simplex_free", "Simplex variable", x_ref);
  Eigen::Index Km1 = x_ref.size() - 1;
  plain_type_t<Vec> y(Km1);
  T stick_len = x_ref.coeff(Km1);
  for (Eigen::Index k = Km1; --k >= 0;) {
    stick_len += x_ref.coeff(k);
    T z_k = x_ref.coeff(k) / stick_len;
    y.coeffRef(k) = logit(z_k) + log(Km1 - k);
    // note: log(Km1 - k) = logit(1.0 / (Km1 + 1 - k));
  }
  return y;
}

/**
 * Overload of `simplex_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto simplex_free(const T& x) {
  return apply_vector_unary<T>::apply(x,
                                      [](auto&& v) { return simplex_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
