#ifndef STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP
#define STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/logit.hpp>
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
 * @tparam T type of elements in the simplex
 * @param x Simplex of dimensionality K.
 * @return Free vector of dimensionality (K-1) that transforms to
 * the simplex.
 * @throw std::domain_error if x is not a valid simplex
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> simplex_free(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::log;

  using size_type = typename index_type<Matrix<T, Dynamic, 1>>::type;

  check_simplex("stan::math::simplex_free", "Simplex variable", x);
  int Km1 = x.size() - 1;
  Eigen::Matrix<T, Eigen::Dynamic, 1> y(Km1);
  T stick_len(x(Km1));
  for (size_type k = Km1; --k >= 0;) {
    stick_len += x(k);
    T z_k(x(k) / stick_len);
    y(k) = logit(z_k) + log(Km1 - k);
    // note: log(Km1 - k) = logit(1.0 / (Km1 + 1 - k));
  }
  return y;
}

}  // namespace math
}  // namespace stan

#endif
