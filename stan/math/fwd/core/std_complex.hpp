#ifndef STAN_MATH_FWD_CORE_STD_COMPLEX_HPP
#define STAN_MATH_FWD_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/meta.hpp>
#include <complex>

namespace std {

/**
 * Specialization of the standard library complex number type for
 * reverse-mode autodiff type `stan::math::fvar<T>`.
 *
 * @tparam T forward-mode autodiff value type
 */
template <typename T>
class complex<stan::math::fvar<T>>
    : public stan::math::complex_base<stan::math::fvar<T>> {
 public:
  using base_t = stan::math::complex_base<stan::math::fvar<T>>;

  /**
   * Construct a complex number with zero real and imaginary parts.
   */
  complex() = default;

  /**
   * Construct a complex number with the specified real part and a zero
   * imaginary part.
   *
   * @tparam Scalar real type (must be assignable to `value_type`)
   * @param[in] re real part
   */
  template <typename U, typename = stan::require_stan_scalar_t<U>>
  complex(const U& re) : base_t(re) {}  // NOLINT(runtime/explicit)

  template <typename U>
  complex(const std::complex<U>& z)  // NOLINT(runtime/explicit)
      : base_t(z.real(), z.imag()) {}

  /**
   * Construct a complex number from the specified real and imaginary
   * parts.
   *
   * @tparam U type of real part
   * @tparam V type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename U, typename V>
  complex(const U& re, const V& im) : base_t(re, im) {}

  /**
   * Set the real and imaginary parts to those of the specified
   * complex number.
   *
   * @tparam U value type of argument
   * @param[in] x complex number to set
   * @return this
   */
  template <typename U, typename = stan::require_arithmetic_t<U>>
  auto& operator=(const std::complex<U>& x) {
    this->re_ = x.real();
    this->im_ = x.imag();
    return *this;
  }
};

}  // namespace std

#endif
