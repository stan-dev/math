#ifndef STAN_MATH_FWD_FUN_LBETA_HPP
#define STAN_MATH_FWD_FUN_LBETA_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>

#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Return digamma(x) - digamma(x + y), the partial derivative of
 * lbeta(x, y) in x. For large x the plain difference loses all digits,
 * so from x = digamma_diff_min_x on it is formed by
 * digamma_diff. Below that value (and for NaN) the plain difference is
 * accurate enough for a gradient and cheaper.
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> lbeta_partial_fwd(const T1& x, const T2& y) {
  if (value_of_rec(x) >= digamma_diff_min_x) {
    return -digamma_diff(x, y);
  }
  return digamma(x) - digamma(x + y);
}
}  // namespace internal

template <typename T>
inline fvar<T> lbeta(const fvar<T>& x1, const fvar<T>& x2) {
  if (value_of_rec(x1.val_) >= internal::digamma_diff_min_x
      || value_of_rec(x2.val_) >= internal::digamma_diff_min_x) {
    return fvar<T>(lbeta(x1.val_, x2.val_),
                   x1.d_ * internal::lbeta_partial_fwd(x1.val_, x2.val_)
                       + x2.d_ * internal::lbeta_partial_fwd(x2.val_, x1.val_));
  }
  // both arguments small: the plain differences share digamma(x1 + x2)
  return fvar<T>(lbeta(x1.val_, x2.val_),
                 x1.d_ * digamma(x1.val_) + x2.d_ * digamma(x2.val_)
                     - (x1.d_ + x2.d_) * digamma(x1.val_ + x2.val_));
}

template <typename T>
inline fvar<T> lbeta(double x1, const fvar<T>& x2) {
  return fvar<T>(lbeta(x1, x2.val_),
                 x2.d_ * internal::lbeta_partial_fwd(x2.val_, x1));
}

template <typename T>
inline fvar<T> lbeta(const fvar<T>& x1, double x2) {
  return fvar<T>(lbeta(x1.val_, x2),
                 x1.d_ * internal::lbeta_partial_fwd(x1.val_, x2));
}
}  // namespace math
}  // namespace stan
#endif
