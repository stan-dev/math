#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

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
   * Constructs complex number from the specified argument, which may
   * be a scalar or complex number.  If the argument is a real type,
   * then it is assigned to the real component and the imaginary
   * component is set to zero; if it is a complex argument, the real
   * and imaginary components are set to their values in the argument.
   *
   * @tparam T type of argument
   * @param[in] re scalar or complex argument
   */
  template <typename T>
  complex(const T& x)  // NOLINT(runtime/explicit)
      : base_t(x) {}

  /**
   * Construct a complex number with zero real part and zero imaginary
   * part.
   */
  complex() : base_t() {}

  /**
   * Construct a complex number with the same real and imaginary
   * components as the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other) : base_t(other) {}

  /**
   * Destroy this complex number.  The implementation is a no-op.
   */
  ~complex() {}

  /**
   * Assign the specified value to this complex number and return this
   * complex number.  If the value is a scalar, assign it to the real
   * component and assign the imaginary component to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number after assigning value
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this and returns this.  If other is a scalar, it
   * adds it to the real component.
   *
   * @tparam V type of value
   * @param[in] x value to add
   * @return this complex number after adding value
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this and returns this.  If other is a
   * scalar, it subtracts it from the real component.
   *
   * @tparam V type of value
   * @param[in] x value to subtract
   * @return this complex number after subtracting value
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Multiplies this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const V& other) {
    return base_t::operator*=(other);
  }

  /**
   * Divides this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const V& other) {
    return base_t::operator/=(other);
  }
};

}  // namespace std

#endif
