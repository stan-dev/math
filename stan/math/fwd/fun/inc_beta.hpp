#ifndef STAN_MATH_FWD_FUN_INC_BETA_HPP
#define STAN_MATH_FWD_FUN_INC_BETA_HPP

#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddb.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/digamma.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, const fvar<T>& b, const fvar<T>& x) {
  T aT = a.val_;
  T bT = b.val_;
  T xT = x.val_;
  T digamma_ab = digamma(aT + bT);
  T d_a = inc_beta_dda(aT, bT, xT, digamma(aT), digamma_ab);
  T d_b = inc_beta_ddb(aT, bT, xT, digamma(bT), digamma_ab);
  T d_x = inc_beta_ddz(aT, bT, xT);
  return fvar<T>(inc_beta(aT, bT, xT), a.d_ * d_a + b.d_ * d_b + x.d_ * d_x);
}

template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T>& b, const fvar<T>& x) {
  T aT = T(a);
  T bT = b.val_;
  T xT = x.val_;
  T digamma_ab = digamma(aT + bT);
  T d_b = inc_beta_ddb(aT, bT, xT, digamma(bT), digamma_ab);
  T d_x = inc_beta_ddz(aT, bT, xT);
  return fvar<T>(inc_beta(aT, bT, xT), b.d_ * d_b + x.d_ * d_x);
}

template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, double b, const fvar<T>& x) {
  T aT = a.val_;
  T bT = T(b);
  T xT = x.val_;
  T digamma_ab = digamma(aT + bT);
  T d_a = inc_beta_dda(aT, bT, xT, digamma(aT), digamma_ab);
  T d_x = inc_beta_ddz(aT, bT, xT);
  return fvar<T>(inc_beta(aT, bT, xT), a.d_ * d_a + x.d_ * d_x);
}

template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, const fvar<T>& b, double x) {
  T aT = a.val_;
  T bT = b.val_;
  T xT = T(x);
  T digamma_ab = digamma(aT + bT);
  T d_a = inc_beta_dda(aT, bT, xT, digamma(aT), digamma_ab);
  T d_b = inc_beta_ddb(aT, bT, xT, digamma(bT), digamma_ab);
  return fvar<T>(inc_beta(aT, bT, xT), a.d_ * d_a + b.d_ * d_b);
}

template <typename T>
inline fvar<T> inc_beta(double a, double b, const fvar<T>& x) {
  T aT = T(a);
  T bT = T(b);
  T xT = x.val_;
  T d_x = inc_beta_ddz(aT, bT, xT);
  return fvar<T>(inc_beta(aT, bT, xT), x.d_ * d_x);
}

template <typename T>
inline fvar<T> inc_beta(const fvar<T>& a, double b, double x) {
  T aT = a.val_;
  T bT = T(b);
  T xT = T(x);
  T digamma_ab = digamma(aT + bT);
  T d_a = inc_beta_dda(aT, bT, xT, digamma(aT), digamma_ab);
  return fvar<T>(inc_beta(aT, bT, xT), a.d_ * d_a);
}

template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T>& b, double x) {
  T aT = T(a);
  T bT = b.val_;
  T xT = T(x);
  T digamma_ab = digamma(aT + bT);
  T d_b = inc_beta_ddb(aT, bT, xT, digamma(bT), digamma_ab);
  return fvar<T>(inc_beta(aT, bT, xT), b.d_ * d_b);
}

}  // namespace math
}  // namespace stan

#endif
