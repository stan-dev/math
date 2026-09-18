#ifndef STAN_MATH_REV_FUN_EIGEN_SPECIALFUNCTIONS_HPP
#define STAN_MATH_REV_FUN_EIGEN_SPECIALFUNCTIONS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <unsupported/Eigen/SpecialFunctions>

namespace Eigen {
namespace internal {

/**
 * Support for Eigen's incomplete gamma routines with Stan's reverse-mode
 * scalar. See `stan/math/fwd/fun/Eigen_SpecialFunctions.hpp` for the
 * reasoning; this is the same three helpers for `var`.
 *
 * No function inside Stan Math reaches `igamma_der_a` with `var`. Every
 * distribution calls the gradient roots with `T_partials_return`, which is
 * `double` for reverse mode. These specializations exist so that a direct
 * call such as `grad_reg_inc_gamma(var, var, var, var)` keeps working after
 * the hand-written series was removed.
 *
 * A `var` call builds one autodiff node per operation inside the series, so
 * it is far more expensive than the `double` call. That was equally true of
 * the previous implementation.
 */
template <>
struct lgamma_impl<stan::math::var> {
  static EIGEN_STRONG_INLINE stan::math::var run(const stan::math::var& x) {
    return stan::math::lgamma(x);
  }
};

template <>
struct digamma_impl<stan::math::var> {
  static EIGEN_STRONG_INLINE stan::math::var run(const stan::math::var& x) {
    return stan::math::digamma(x);
  }
};

/**
 * The Cephes constants are properties of the underlying floating point
 * format, so they are constants with no adjoint.
 */
template <>
struct cephes_helper<stan::math::var> {
  static EIGEN_STRONG_INLINE stan::math::var machep() {
    return stan::math::var(cephes_helper<double>::machep());
  }
  static EIGEN_STRONG_INLINE stan::math::var big() {
    return stan::math::var(cephes_helper<double>::big());
  }
  static EIGEN_STRONG_INLINE stan::math::var biginv() {
    return stan::math::var(cephes_helper<double>::biginv());
  }
};

}  // namespace internal
}  // namespace Eigen

#endif
