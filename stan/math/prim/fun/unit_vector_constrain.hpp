#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector x.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam Vec type with a defined `operator[]`
 * @param x vector of K unrestricted variables
 * @return Unit length vector of dimension K
 */
template <typename Vec, require_vector_like_t<Vec>* = nullptr>
auto unit_vector_constrain(Vec&& x) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "x", x);
  value_type_t<Vec> SN = dot_self(x);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  return divide(x, sqrt(SN));
}

/**
 * Return the unit length vector corresponding to the free vector x.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 *
 * @tparam Vec type with a defined `operator[]`
 * @tparam T type of log probability.
 * @param x vector of K unrestricted variables
 * @param lp Log probability reference to increment.
 * @return Unit length vector of dimension K
 */
template <typename Vec, typename T, require_vector_like_t<Vec>* = nullptr>
auto unit_vector_constrain(Vec&& x, T& lp) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "x", x);
  value_type_t<Vec> SN = dot_self(x);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  lp -= 0.5 * SN;
  return divide(x, sqrt(SN));
}

}  // namespace math
}  // namespace stan

#endif
