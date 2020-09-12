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
 * @tparam ColVec type of the simplex (must be a column vector)
 * @param x Simplex of dimensionality K.
 * @return Free vector of dimensionality (K-1) that transforms to
 * the simplex.
 * @throw std::domain_error if x is not a valid simplex
 */
template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
auto simplex_free(const ColVec& x) {
  using std::log;
  using T = value_type_t<ColVec>;

  check_simplex("stan::math::simplex_free", "Simplex variable", x);
  int Km1 = x.size() - 1;
  Eigen::Matrix<T, Eigen::Dynamic, 1> y(Km1);
  T stick_len = x.coeff(Km1);
  for (Eigen::Index k = Km1; --k >= 0;) {
    stick_len += x.coeff(k);
    T z_k = x.coeff(k) / stick_len;
    y.coeffRef(k) = logit(z_k) + log(Km1 - k);
    // note: log(Km1 - k) = logit(1.0 / (Km1 + 1 - k));
  }
  return y;
}

}  // namespace math
}  // namespace stan

#endif
