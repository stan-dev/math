#ifndef STAN_MATH_FWD_FUN_INC_BETA_HPP
#define STAN_MATH_FWD_FUN_INC_BETA_HPP

#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/beta.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/pow.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, const fvar<T>& b, const fvar<T>& x) {
  using std::pow;
  T d_a;
  T d_b;
  const T beta_ab = beta(a.val_, b.val_);
  grad_reg_inc_beta(d_a, d_b, a.val_, b.val_, x.val_, digamma(a.val_),
                    digamma(b.val_), digamma(a.val_ + b.val_), beta_ab);
  T d_x = pow((1 - x.val_), b.val_ - 1) * pow(x.val_, a.val_ - 1) / beta_ab;
  return fvar<T>(inc_beta(a.val_, b.val_, x.val_),
                 a.d_ * d_a + b.d_ * d_b + x.d_ * d_x);
}

template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T>& b, const fvar<T>& x) {
  return inc_beta(fvar<T>(a), b, x);
}
template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, double b, const fvar<T>& x) {
  return inc_beta(a, fvar<T>(b), x);
}
template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, const fvar<T>& b, double x) {
  return inc_beta(a, b, fvar<T>(x));
}

template <typename T>
inline fvar<T> inc_beta(double a, double b, const fvar<T>& x) {
  return inc_beta(fvar<T>(a), fvar<T>(b), x);
}
template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, double b, double x) {
  return inc_beta(a, fvar<T>(b), fvar<T>(x));
}
template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T>& b, double x) {
  return inc_beta(fvar<T>(a), b, fvar<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
