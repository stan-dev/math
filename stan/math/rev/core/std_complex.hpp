#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/rev/core/var.hpp>
#include <cmath>
#include <complex>

namespace std {

/**
 * Specialization of the standard library complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class complex<stan::math::var>
    : public stan::math::complex_base<stan::math::var> {
 public:
  using base_t = stan::math::complex_base<stan::math::var>;

  /**
   * Construct a complex number with zero real and imaginary parts.
   */
  complex() = default;

  /**
   * Construct a complex number from real and imaginary parts.
   *
   * @tparam U type of real part (assignable to `value_type`)
   * @tparam V type of imaginary part (assignable to `value_type`)
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename U, typename V>
  complex(const U& re, const V& im) : base_t(re, im) {}

  /**
   * Construct a complex number with specified real part and zero
   * imaginary part.
   *
   * @tparam U type of real part (assignable to `value_type`)
   * @param[in] re real part
   */
  template <typename U, typename = stan::require_stan_scalar_t<U>>
  complex(const U& re) : base_t(re) {}  // NOLINT(runtime/explicit)

  template <typename U>
  complex(const std::complex<U>& z)  // NOLINT(runtime/explicit)
      : base_t(z.real(), z.imag()) {}

  /**
   * Set the real and imaginary components of this complex number to
   * those of the specified complex number.
   *
   * @tparam U value type of argument (must be an arithmetic type)
   * @param[in] x complex argument
   * @return this
   */
  template <typename U, typename = stan::require_arithmetic_t<U>>
  auto& operator=(const std::complex<U>& x) {
    re_ = x.real();
    im_ = x.imag();
    return *this;
  }
};

}  // namespace std

#endif
