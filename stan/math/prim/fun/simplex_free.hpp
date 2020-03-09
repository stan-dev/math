#ifndef STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP
#define STAN_MATH_PRIM_FUN_SIMPLEX_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
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
 * @tparam Vec type deriving from `Eigen::MatrixBase` with rows or columns
 * equal to 1.
 * @param x Simplex of dimensionality K.
 * @return Free vector of dimensionality (K-1) that transforms to
 * the simplex.
 * @throw std::domain_error if x is not a valid simplex
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
auto simplex_free(Vec&& x) {
  using std::log;
  check_simplex("stan::math::simplex_free", "Simplex variable", x);
  auto Km1 = x.size() - 1;
  auto ret_size = Km1 >= 0 ? Km1 : 0;
  plain_type_t<Vec> y(ret_size);
  if (y.size() == 0) {
    return y;
  }
  auto stick_len = x[Km1];
  for (auto k = Km1 - 1; k > -1; k--) {
    stick_len += x[k];
    auto z_k = x[k] / stick_len;
    y[k] = logit(z_k) + log(Km1 - k);
    // note: log(Km1 - k) = logit(1.0 / (Km1 + 1 - k));
  }
  return y;
}

}  // namespace math
}  // namespace stan

#endif
