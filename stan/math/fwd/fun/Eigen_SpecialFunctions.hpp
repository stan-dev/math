#ifndef STAN_MATH_FWD_FUN_EIGEN_SPECIALFUNCTIONS_HPP
#define STAN_MATH_FWD_FUN_EIGEN_SPECIALFUNCTIONS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <unsupported/Eigen/SpecialFunctions>

namespace Eigen {
namespace internal {

/**
 * Support for Eigen's incomplete gamma routines with Stan's forward-mode
 * scalars.
 *
 * `Eigen::internal::igamma_generic_impl` is written for a general scalar
 * type and propagates the derivative through the Cephes power series and
 * the Cephes continued fraction. Only three of its helpers are restricted
 * to `float` and `double`. Supplying those three here lets
 * `Eigen::numext::igamma_der_a` serve `fvar<T>` at every autodiff order,
 * so `grad_reg_inc_gamma` and `grad_reg_lower_inc_gamma` can use one
 * algorithm for every scalar type.
 *
 * This depends on the names `lgamma_impl`, `digamma_impl` and
 * `cephes_helper` in `Eigen::internal`. An Eigen upgrade that renames them
 * fails at compile time, not silently.
 */
template <typename T>
struct lgamma_impl<stan::math::fvar<T>> {
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE stan::math::fvar<T> run(
      const stan::math::fvar<T>& x) {
    return stan::math::lgamma(x);
  }
};

template <typename T>
struct digamma_impl<stan::math::fvar<T>> {
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE stan::math::fvar<T> run(
      const stan::math::fvar<T>& x) {
    return stan::math::digamma(x);
  }
};

/**
 * The Cephes constants are properties of the underlying floating point
 * format, not of the autodiff type, so they carry a zero tangent.
 */
template <typename T>
struct cephes_helper<stan::math::fvar<T>> {
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE stan::math::fvar<T> machep() {
    return stan::math::fvar<T>(cephes_helper<double>::machep());
  }
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE stan::math::fvar<T> big() {
    return stan::math::fvar<T>(cephes_helper<double>::big());
  }
  EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE stan::math::fvar<T> biginv() {
    return stan::math::fvar<T>(cephes_helper<double>::biginv());
  }
};

}  // namespace internal
}  // namespace Eigen

#endif
