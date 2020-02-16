#ifndef STAN_MATH_PRIM_CORE_COMPLEX_BASE_HPP
#define STAN_MATH_PRIM_CORE_COMPLEX_BASE_HPP

#include <stan/math/prim/fun/square.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Base class for complex numbers.  The template parameter `V` is the
 * value type.  This class is used like the curiously recursive
 * template pattern (CRTP) with derived type `complex<V>`.  Thus any
 * extending class must be of the form `complex<V>`.
 *
 * @tparam V value type for extending complex class
 */
template <typename V>
class complex_base {
 protected:
  /**
   * Real part.
   */
  V re_;

  /**
   * Imaginary part.
   */
  V im_;

 public:
  /**
   * Type of real and imaginary parts.
   */
  using value_type = V;

  /**
   * Derived complex type used for function return types.
   */
  using complex_type = std::complex<value_type>;

  /**
   * Construct a complex base with the specified real and complex
   * parts.
   *
   * @tparam T real type (must be assignable to `V`)
   * @tparam U imaginary type (must be assignable to `V`)
   * @param[in] re real part
   * @param[in] im imaginary part (default 0)
   */
  template <typename T, typename U>
  complex_base(const T& re, const U& im = U(0))  // NOLINT(runtime/explicit)
      : re_(re), im_(im) {}

  /**
   * Constructs complex number from real and imaginary parts.
   *
   * @param[in] re the real part
   * @param[in] im the imaginary part
   */
  complex_base(const value_type& re = value_type(0),
               const value_type& im = value_type(0))
      : re_(re), im_(im) {}

  /**
   * Constructs the complex number from the specified complex number
   * of a different type.
   *
   * @param[in] other another complex to use as source
   */
  template <typename T>
  complex_base(const std::complex<T>& other)  // NOLINT(runtime/explicit)
      : complex_base(other.real(), other.imag()) {}

  /**
   * Destroy this complex number.
   */
  ~complex_base() {}

  /**
   * Return a reference to thi class cast to the derived complex
   * type.
   *
   * @return reference to this class cast to the complex return type
   */
  complex_type& derived_complex_ref() {
    return static_cast<complex_type&>(*this);
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  complex_type& operator=(const V& x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  complex_type& operator=(float x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  complex_type& operator=(double x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  complex_type& operator=(long double x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @tparam T type of value, which must be assignable to the complex
   * value type
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename T,
            typename
            = typename std::enable_if_t<std::is_assignable<V, T>::value>>
  complex_type& operator=(const T& x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam T value type of argument (must be assignable to
   * `value_type`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename T>
  complex_type& operator=(const std::complex<T>& x) {
    re_ = x.real();
    im_ = x.imag();
    return derived_complex_ref();
  }

  complex_type& operator=(const std::complex<V>& x) {
    re_ = x.real();
    im_ = x.imag();
    return derived_complex_ref();
  }

  /**
   * Return the real part.
   *
   * @return the real part
   */
  value_type real() const { return re_; }

  /**
   * Set the real part to the specified value.
   *
   * @param[in] x the value to set the real part to
   */
  void real(const value_type& x) { re_ = x; }

  /**
   * Return the imaginary part.
   *
   * @return the imaginary part
   */
  value_type imag() const { return im_; }

  /**
   * Set the imaginary part to the specified value.
   *
   * @param[in] x the value to set the imaginary part to
   */
  void imag(const value_type& x) { im_ = x; }

  /**
   * Adds other to this.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const X& other) {
    re_ += other;
    return derived_complex_ref();
  }

  /**
   * Adds other to this.
   *
   * @tparam value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const std::complex<X>& other) {
    re_ += other.real();
    im_ += other.imag();
    return derived_complex_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X type of scalar argument (assignable to vlaue type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const X& other) {
    re_ -= other;
    return derived_complex_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const std::complex<X>& other) {
    re_ -= other.real();
    im_ -= other.imag();
    return derived_complex_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const X& other) {
    re_ *= other;
    im_ *= other;
    return derived_complex_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const std::complex<X>& other) {
    value_type re_temp = re_ * other.real() - im_ * other.imag();
    im_ = re_ * other.imag() + other.real() * im_;
    re_ = re_temp;
    return derived_complex_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator/=(const X& other) {
    re_ /= other;
    im_ /= other;
    return derived_complex_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator/=(const std::complex<X>& other) {
    using stan::math::square;
    value_type sum_sq_im = square(other.real()) + square(other.imag());
    value_type re_temp = (re_ * other.real() + im_ * other.imag()) / sum_sq_im;
    im_ = (im_ * other.real() - re_ * other.imag()) / sum_sq_im;
    re_ = re_temp;
    return derived_complex_ref();
  }
};

}  // namespace math
}  // namespace stan

#endif
