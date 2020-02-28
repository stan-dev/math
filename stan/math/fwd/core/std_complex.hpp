#ifndef STAN_MATH_FWD_CORE_STD_COMPLEX_HPP
#define STAN_MATH_FWD_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <complex>

namespace std {

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::fvar<T>`.
 */
template <typename T>
class complex<stan::math::fvar<T>>
    : public stan::math::complex_base<stan::math::fvar<T>> {
 public:
  using base_t = stan::math::complex_base<stan::math::fvar<T>>;
  using value_type = stan::math::fvar<T>;
  using complex_type = complex<value_type>;

  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam U type of real part
   * @tparam V type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename U, typename V>
  complex(const U& re, const V& im) : base_t(re, im) {}

  /**
   * Constructs complex number from argument, which may be a scalar or
   * a complex number.  If it is a scalar, the real part is set to the
   * scalar value and the imaginary part to zero;  if it is complex,
   * the real and imaginary parts are set to those of the argument.
   *
   * @tparam U type of argument
   * @param[in] x the scalar or complex argument
   */
  template <typename U>
  complex(const U& x) : base_t(x) {}  // NOLINT(runtime/explicit)

  /**
   * Construct a complex number with zero real and imaginary
   * components.
   */
  complex() : base_t() {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex
   * number and set imaginary part to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number after assigning other
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this and return this.
   *
   * @tparam V type of value
   * @param[in] other value to add
   * @return this complex number after adding argument
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this and return this.
   *
   * @tparam V type of value
   * @param[in] other value to subtract
   * @return this complex number after subtracting argument
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Multiplies this by other and return this.
   *
   * @tparam V type of value
   * @param[in] other value to multiply
   * @return this complex number after multiplying by argument
   */
  template <typename V>
  complex_type& operator*=(const V& other) {
    return base_t::operator*=(other);
  }

  /**
   * Divides this by other and return this.
   *
   * @tparam V type of value
   * @param[in] other value to divide by
   * @return this complex number after dividing by argument
   */
  template <typename V>
  complex_type& operator/=(const V& other) {
    return base_t::operator/=(other);
  }
};

}  // namespace std

#endif
