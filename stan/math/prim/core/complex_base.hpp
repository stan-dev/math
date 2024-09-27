#ifndef STAN_MATH_PRIM_CORE_COMPLEX_HPP
#define STAN_MATH_PRIM_CORE_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Base class for complex numbers.
 *
 * @tparam ValueType type of real and imaginary parts
 */
template <typename ValueType>
class complex {
 public:
  /**
   * Type of real and imaginary parts
   */
  using value_type = ValueType;
  using complex_type = complex<value_type>;

  /**
   * Construct a complex base with zero real and imaginary parts.
   */
  complex() = default;

  /**
   * Construct a complex base with the specified real part and a zero
   * imaginary part.
   *
   * @tparam U real type (assignable to `value_type`)
   * @param[in] re real part
   */
  template <typename U>              // , typename = require_stan_scalar_t<U>>
  complex(const U& re) : re_(re) {}  // NOLINT(runtime/explicit)

  /**
   * Construct a complex base with the specified real and imaginary
   * parts.
   *
   * @tparam U real type (assignable to `value_type`)
   * @tparam V imaginary type (assignable to `value_type`)
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename U, typename V>
  complex(const U& re, const V& im) : re_(re), im_(im) {}

  /**
   * Return the real part.
   *
   * @return real part
   */
  value_type real() const { return re_; }

  /**
   * Set the real part to the specified value.
   *
   * @param[in] re real part
   */
  void real(const value_type& re) { re_ = re; }

  /**
   * Return the imaginary part.
   *
   * @return imaginary part
   */
  value_type imag() const { return im_; }

  /**
   * Set the imaginary part to the specified value.
   *
   * @param[in] im imaginary part
   */
  void imag(const value_type& im) { im_ = im; }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @tparam U argument type (assignable to `value_type`)
   * @param[in] re real part
   * @return this
   */
  template <typename U, typename = require_stan_scalar_t<U>>
  complex_type& operator=(U&& re) {
    re_ = re;
    im_ = 0;
    return *this;
  }

  /**
   * Add specified real value to real part.
   *
   * @tparam U argument type (assignable to `value_type`)
   * @param[in] x real number to add
   * @return this
   */
  template <typename U>
  complex_type& operator+=(const U& x) {
    re_ += x;
    return *this;
  }

  /**
   * Adds specified complex number to this.
   *
   * @tparam U value type of argument (assignable to `value_type`)
   * @param[in] other complex number to add
   * @return this
   */
  template <typename U>
  complex_type& operator+=(const complex<U>& other) {
    re_ += other.real();
    im_ += other.imag();
    return *this;
  }

  /**
   * Subtracts specified real number from real part.
   *
   * @tparam U argument type (assignable to `value_type`)
   * @param[in] x real number to subtract
   * @return this
   */
  template <typename U>
  complex_type& operator-=(const U& x) {
    re_ -= x;
    return *this;
  }

  /**
   * Subtracts specified complex number from this.
   *
   * @tparam U value type of argument (assignable to `value_type`)
   * @param[in] other complex number to subtract
   * @return this
   */
  template <typename U>
  complex_type& operator-=(const complex<U>& other) {
    re_ -= other.real();
    im_ -= other.imag();
    return *this;
  }

  /**
   * Multiplies this by the specified real number.
   *
   * @tparam U type of argument (assignable to `value_type`)
   * @param[in] x real number to multiply
   * @return this
   */
  template <typename U>
  complex_type& operator*=(const U& x) {
    re_ *= x;
    im_ *= x;
    return *this;
  }

  /**
   * Multiplies this by specified complex number.
   *
   * @tparam U value type of argument (assignable to `value_type`)
   * @param[in] other complex number to multiply
   * @return this
   */
  template <typename U>
  complex_type& operator*=(const complex<U>& other) {
    value_type re_temp = re_ * other.real() - im_ * other.imag();
    im_ = re_ * other.imag() + other.real() * im_;
    re_ = re_temp;
    return *this;
  }

  /**
   * Divides this by the specified real number.
   *
   * @tparam U type of argument (assignable to `value_type`)
   * @param[in] x real number to divide by
   * @return this
   */
  template <typename U>
  complex_type& operator/=(const U& x) {
    re_ /= x;
    im_ /= x;
    return *this;
  }

  /**
   * Divides this by the specified complex number.
   *
   * @tparam U value type of argument (assignable to `value_type`)
   * @param[in] other number to divide by
   * @return this
   */
  template <typename U>
  complex_type& operator/=(const complex<U>& other) {
    value_type sum_sq_im
        = (other.real() * other.real()) + (other.imag() * other.imag());
    value_type re_temp = (re_ * other.real() + im_ * other.imag()) / sum_sq_im;
    im_ = (im_ * other.real() - re_ * other.imag()) / sum_sq_im;
    re_ = re_temp;
    return *this;
  }

  template <typename V = ValueType, require_arithmetic_t<V>* = nullptr>
  stan::math::complex<V> std_complex() const {
    return static_cast<stan::math::complex<V>>(*this);
  }

 protected:
  /**
   * Real part
   */
  value_type re_{0};

  /**
   * Imaginary part
   */
  value_type im_{0};
};

}  // namespace math
}  // namespace stan

#endif
