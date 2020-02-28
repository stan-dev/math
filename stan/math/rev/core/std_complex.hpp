#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/rev/core/var.hpp>
#include <cmath>
#include <complex>

namespace std {

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class complex<stan::math::var>
    : public stan::math::complex_base<stan::math::var> {
 public:
  using base_t = stan::math::complex_base<stan::math::var>;
  using value_type = stan::math::var;
  using complex_type = complex<value_type>;
  complex() = default;
  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam T type of real part
   * @tparam U type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename T, typename U>
  complex(const T& re, const U& im) : base_t(re, im) {}

  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam T type of real part
   * @tparam U type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename T>
  complex(const T& re) : base_t(re) {}


};

}  // namespace std

#endif
