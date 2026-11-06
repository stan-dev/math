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
 * The same three `Eigen::internal` helper specializations as
 * `stan/math/fwd/fun/Eigen_SpecialFunctions.hpp`, for `var`, so that
 * `Eigen::numext::igamma_der_a` accepts `var` arguments.
 *
 * Inside Stan Math the gradient roots are called with `double` partials in
 * reverse mode; this path serves direct calls with `var`. Such a call
 * records one autodiff node per operation of the series, so it is far more
 * expensive than the `double` call.
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
