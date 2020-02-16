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
  using base_t = complex_base<stan::math::var>;

  /**
   * Constructs complex number from real part or with default zero
   * value, setting imaginary part to zero.  This is the nullary
   * constructor when the default of zero is used for the real
   * component.
   *
   * @param[in] re the real part (defaults to zero)
   */
  complex(const value_type& re = value_type(0))  // NOLINT(runtime/explicit)
      : complex_base(re) {}

  /**
   * Construct a complex number with the same real and imaginary
   * components as the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other) : complex_base(other) {}

  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam V1 type of real part
   * @tparam V2 type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename V1, typename V2>
  complex(const V1& re, const V2& im) : complex_base(re, im) {}

  /**
   * Construct a complex number with the specified real component and
   * zero imaginary component.
   *
   * @tparam T arithmetic type of argument
   * @param[in] x real component
   */
  template <typename T,
            typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  complex(T x) : complex_base(x) {}  // NOLINT(runtime/explicit)

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
