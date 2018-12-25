#ifndef STAN_MATH_PRIM_MAT_FUN_OFFSET_MULTIPLIER_FREE_HPP
#define STAN_MATH_PRIM_MAT_FUN_OFFSET_MULTIPLIER_FREE_HPP

#include <stan/math/prim/mat/err/check_cholesky_factor.hpp>
#include <stan/math/prim/mat/err/check_finite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <boost/math/tools/promotion.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the unconstrained vector that transforms to the
 * specified offset and multiplier constrained vector given the specified
 * offset and multiplier.
 *
 * <p>The transfrom in <code>offset_multiplier_constrain(T, double,
 * double)</code>, is reversed by the reverse affine transformation,
 *
 * <p>\f$f^{-1}(y) = \frac{y - L}{S}\f$
 *
 * where \f$L\f$ and \f$S\f$ are the offset and multiplier.
 *
 * <p>If the offset is zero and multiplier is one,
 * this function reduces to  <code>identity_free(y)</code>.
 *
 * @tparam T type of vector
 * @tparam L type of offset
 * @tparam S type of multiplier
 * @param y constrained value
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return the free vector that transforms to the input vector
 *   given the offset and multiplier
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename L, typename S>
inline Eigen::Matrix<typename return_type<T, L, S>::type, -1, 1>
offset_multiplier_free(const Eigen::Matrix<T, -1, 1>& y,
                       const Eigen::Matrix<L, -1, 1>& mu,
                       const Eigen::Matrix<S, -1, -1>& sigma) {
  static const char* function = "offset_multiplier_free";
  check_finite(function, "offset", mu);
  check_finite(function, "multiplier", sigma);
  check_cholesky_factor(function, "multiplier", sigma);
  check_square(function, "multiplier", sigma);
  check_consistent_sizes(function, "multiplier", sigma.col(0),
                         "contrained vector", y);
  const size_t N = sigma.col(0).size();
  if (sigma == Eigen::Matrix<S, -1, -1>::Identity(N, N)) {
    if (mu == Eigen::Matrix<L, -1, 1>::Zero(N, 1))
      return y;
    return y - mu;
  }
  return mdivide_left_ldlt(sigma, y - mu);
}

}  // namespace math
}  // namespace stan
#endif
