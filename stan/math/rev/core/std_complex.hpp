#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/scal/fun/atan2.hpp>
#include <stan/math/rev/scal/fun/cos.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/rev/scal/fun/hypot.hpp>
#include <stan/math/rev/scal/fun/is_inf.hpp>
#include <stan/math/rev/scal/fun/is_nan.hpp>
#include <stan/math/rev/scal/fun/log.hpp>
#include <stan/math/rev/scal/fun/sin.hpp>
#include <stan/math/rev/scal/fun/square.hpp>
#include <stan/math/rev/scal/fun/sqrt.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <limits>

#include <iostream>

// std SPECIALIZATIONS FOR STAN W/O COMPLEX
// ===================================================================

namespace std {
/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <>
struct iterator_traits<stan::math::var> {
  /**
   * Iterator category for traits.
   */
  typedef random_access_iterator_tag iterator_category;

  /**
   * Type for difference between pointers.
   */
  typedef ptrdiff_t difference_type;

  /**
   * Type for value of pointer to values.
   */
  typedef stan::math::var value_type;

  /**
   * Type of pointer to variables.
   */
  typedef stan::math::var* pointer;

  /**
   * Type of reference to variables.
   */
  typedef stan::math::var& reference;
};
}  // namespace std

// STAN CODE W/O complex<var>
// ===================================================================

namespace stan {
namespace math {
/**
 * Return `true` if the specified argument is negative and `false`
 * otherwise.
 *
 * Overloads `std::signbit` from `<cmath>`.
 *
 * @tparam T type of argument
 * @param v argument
 * @return `true` if the argument is negative
 */
template <typename T>
bool signbit(const T& v) {
  using std::signbit;
  return signbit(v.val());
}

/**
 * Return true if specified argument is infinite.
 *
 * Overloads `std::isinf` from `<cmath>`.
 *
 * @tparam T type of argument
 * @param v argument
 * @return true if argument is infinite
 */
template <typename T>
bool isinf(const T& v) {
  using std::isinf;
  return isinf(v.val());
}

/**
 * Return true if specified argument is finite.
 *
 * Overloads `std::isfinite` from `<cmath>`.
 *
 * @tparam T type of argument
 * @param v argument
 * @return true if argument is finite
 */
template <typename T>
bool isfinite(const T& v) {
  using std::isfinite;
  return isfinite(v.val());
}

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>`.
 *
 * @tparam T type of argument
 * @param v argument
 * @return true if argument is not-a-number
 */
template <typename T>
bool isnan(const T& v) {
  using std::isnan;
  return isnan(v.val());
}

/**
 * Return true if specified argument is normal.  A number is normal if
 * it is finite, non-zero and not subnormal.
 *
 * Overloads `std::isnormal` from `<cmath>`.
 *
 * @tparam T type of argument
 * @param v argument
 * @return true if argument is normal
 */
template <typename T>
bool isnormal(const T& v) {
  using std::isnormal;
  return isnormal(v.val());
}

/**
 * Return the negation of the first argument if the first and second
 * argument have different signs, otherwise return a copy of the first
 * argument.
 *
 * Overload of `std::copysign` from `cmath`.
 *
 * @tparam T type of first argument
 * @tparam U type of second argument
 * @param x first complex argument
 * @param y second complex argument
 * @return copy of second argument, negated if necessary to match sign
 * of first argument
 */
template <typename T, typename U>
T copysign(const T& x, const U& y) {
  // 0 is considered positive
  return (x < 0 && y >= 0) || (x > 0 && y < 0) ? -x : x;
}

/**
 * Return the negation of the first argument if the first and second
 * argument have different signs, otherwise return a copy of the first
 * argument.
 *
 * This one is not an overload of a standard library function; it's
 * for utility in the complex overloads for `var` types.
 *
 * @tparam T value type of first argument
 * @tparam U value type of second argument
 * @param x first complex argument
 * @param y second complex argument
 * @return copy of second argument, with components negated if
 * necessary to match sign of first argument
 */
template <typename T, typename U>
std::complex<T> copysign(const std::complex<T>& y, const std::complex<U>& x) {
  return {copysign(y.real(), x.real()), copysign(y.imag(), x.imag())};
}

/**
 * Return the specified complex number multiplied by `i`.
 *
 * @tparam value type of complex argument
 * @param z complex argument
 * @return argument multipled by `i`
 */
template <typename T>
std::complex<T> i_times(const std::complex<T>& z) {
  return {-z.imag(), z.real()};
}

/**
 * Return the specified complex number multiplied by `-i`.
 *
 * @tparam value type of complex argument
 * @param z complex argument
 * @return argument multipled by `-i`
 */
template <typename T>
std::complex<T> neg_i_times(const std::complex<T>& z) {
  return {z.imag(), -z.real()};
}

/**
 * Return the complex number with values the same as the specified
 * complex argument, but one lower type.  For example, an argument
 * type of `complex<var>` is converted to a return type
 * `complex<double>` and `complex<fvar<fvar<var>>` to
 * `complex<fvar<var>>`.
 *
 * @tparam T value type of complex argument
 * @param[in] complex argument
 * @return complex number with same value as argument
 */
template <typename T>
std::complex<typename T::Scalar> value_of(const std::complex<T>& z) {
  return {z.real().val(), z.imag().val()};
}
}  // namespace math
}  // namespace stan

namespace stan {
namespace math {
/**
 * CRTP base calss for complex numbers.  Rather than typical CRTP, the
 * template variable is just the value type for the complex number.
 * The extending class will then be `complex<V>`.
 *
 * @tparam V value type for extending complex class
 */
template <typename V>
class complex_base {
 protected:
  V re_;
  V im_;

 public:
  /**
   * Type of real and imaginary parts.
   */
  typedef V value_type;

  typedef std::complex<value_type> complex_type;

  template <typename T, typename U>
  complex_base(const T& x, const U& y = U(0))
      : complex_base(value_type{x}, value_type{y}) {}

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
  complex_base(const std::complex<T>& other)
      : complex_base(other.real(), other.imag()) {}

  /**
   * Destroy this complex number.
   */
  ~complex_base() {}

  complex_type& return_ref() { return static_cast<complex_type&>(*this); }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename T>
  complex_type& operator=(const T& x) {
    re_ = x;
    im_ = 0;
    return return_ref();
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
    return return_ref();
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
   * @tparam X type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const X& other) {
    re_ += other;
    return return_ref();
  }

  /**
   * Adds other to this.
   *
   * @tparam value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const std::complex<X>& other) {
    re_ += other.real();
    im_ += other.imag();
    return return_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const X& other) {
    re_ -= other;
    return return_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const std::complex<X>& other) {
    re_ -= other.real();
    im_ -= other.imag();
    return return_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const X& other) {
    re_ *= other;
    im_ *= other;
    return return_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const std::complex<X>& other) {
    value_type re_temp = re_ * other.real() - im_ * other.imag();
    im_ = re_ * other.imag() + other.real() * im_;
    re_ = re_temp;
    return return_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator/=(const X& other) {
    re_ /= other;
    im_ /= other;
    return return_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X value type of complex argument (assignable to
   * `stan::math::var`)
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
    return return_ref();
  }
};

}  // namespace math
}  // namespace stan

// SPECIALIZATION complex<var> for std::complex
namespace std {

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class complex<stan::math::var>
    : public stan::math::complex_base<stan::math::var> {
 public:
  typedef complex_base<stan::math::var> base_t;

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
   * Constructs complex number from real part or with default zero
   * value, setting imaginary part to zero.
   *
   * @param[in] re the real part
   */
  complex(const value_type& re = value_type(0)) : complex_base(re) {}

  /**
   * Constructs the complex number from the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other) : complex_base(other) {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam V value type of argument (assignable to `stan::math::var`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const complex<V>& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Adds other to this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const complex<V>& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const complex<V>& other) {
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
   * Multiplies this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const complex<V>& other) {
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

  /**
   * Divides this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const complex<V>& other) {
    return base_t::operator/=(other);
  }
};

// After here, it's specializations of function templates; for info, see:
// Walter E. Brown.  2017.  Thou Shalt Not Specialize std Function Templates!
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2017/p0551r1.pdf

/**
 * Return the value of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return a copy of the argument
 */
template <>
complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& val) {
  return val;
}

/**
 * Return the negation of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return negated argument
 */
template <>
complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& val) {
  return complex<stan::math::var>(-val.real(), -val.imag());
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs scalar second argument
 * @return sum of the arguments
 */
template <>
complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}

/**
 * Return the sum of the two arguments.
 *
 * @param[in] lhs scalar first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
complex<stan::math::var> operator+<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs complex second argument
 * @return difference between first and second argument
 */
template <>
complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs scalar second argument
 * @return difference between first and second argument
 */
template <>
complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs scalar first argument
 * @param[in] rhs complex second argument
 * @return difference between first and second argument
 */
template <>
complex<stan::math::var> operator-<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return product of the first and second argument
 */
template <>
complex<stan::math::var> operator*<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return product of the first and second argument
 */
template <>
complex<stan::math::var> operator*<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return product of the first and second argument
 */
template <>
complex<stan::math::var> operator*<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return quotient of the first and second argument
 */
template <>
complex<stan::math::var> operator/<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return quotient of the first and second argument
 */
template <>
complex<stan::math::var> operator/<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return quotient of the first and second argument
 */
template <>
complex<stan::math::var> operator/(const stan::math::var& lhs,
                                   const complex<stan::math::var>& rhs) {
  complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are equal
 */
template <>
bool operator==<stan::math::var>(const complex<stan::math::var>& lhs,
                                 const complex<stan::math::var>& rhs) {
  return lhs.real() == rhs.real() && lhs.imag() == rhs.imag();
}

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return if the arguments are equal
 */
template <>
bool operator==<stan::math::var>(const complex<stan::math::var>& lhs,
                                 const stan::math::var& rhs) {
  return lhs.real() == rhs && lhs.imag() == 0;
}

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are equal
 */
template <>
bool operator==<stan::math::var>(const stan::math::var& lhs,
                                 const complex<stan::math::var>& rhs) {
  return lhs == rhs.real() && 0 == rhs.imag();
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are not equal
 */
template <>
bool operator!=<stan::math::var>(const complex<stan::math::var>& lhs,
                                 const complex<stan::math::var>& rhs) {
  return !(lhs.real() == rhs.real() && lhs.imag() == rhs.imag());
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return if the arguments are not equal
 */
template <>
bool operator!=<stan::math::var>(const complex<stan::math::var>& lhs,
                                 const stan::math::var& rhs) {
  return !(lhs.real() == rhs && lhs.imag() == 0);
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are not equal
 */
template <>
bool operator!=<stan::math::var>(const stan::math::var& lhs,
                                 const complex<stan::math::var>& rhs) {
  return !(lhs == rhs.real() && 0 == rhs.imag());
}

/**
 * Return the real component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
stan::math::var real<stan::math::var>(const complex<stan::math::var>& z) {
  return z.real();
}

/**
 * Return the imaginary component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
stan::math::var imag<stan::math::var>(const complex<stan::math::var>& z) {
  return z.imag();
}

/**
 * Return the magnitude of the specified complex number.
 *
 * @param[in] complex argument
 * @return magnitude of argument
 */
template <>
stan::math::var abs<stan::math::var>(const complex<stan::math::var>& z) {
  return hypot(real(z), imag(z));
}

/**
 * Return the phase angle of the specified complex number.
 *
 * @param[in] z complex argument
 * @return phase angle of argument
 */
template <>
stan::math::var arg<stan::math::var>(const complex<stan::math::var>& z) {
  return atan2(imag(z), real(z));
}

/**
 * Return the squared magnitude of the specified complex number.
 *
 * @param[in] z complex argument
 * @return squared magnitude of argument
 */
template <>
stan::math::var norm<stan::math::var>(const complex<stan::math::var>& z) {
  return square(real(z)) + square(imag(z));
}

/**
 * Return the complex conjugate of the specified complex number.
 *
 * @param[in] z complex argument
 * @return complex conjugate of argument
 */
template <>
complex<stan::math::var> conj<stan::math::var>(
    const complex<stan::math::var>& z) {
  return {z.real(), -z.imag()};
}

/**
 * Return the projection of the specified complex arguent onto the
 * Riemann sphere.
 *
 * @param[in] z complex argument
 * @return project of argument onto the Riemann sphere
 */
template <>
complex<stan::math::var> proj<stan::math::var>(
    const complex<stan::math::var>& z) {
  if (stan::math::is_inf(z.real()) || stan::math::is_inf(z.imag())) {
    return {numeric_limits<stan::math::var>::infinity(),
            z.imag() < 0 ? -0.0 : 0.0};
  }
  return z;
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <>
complex<stan::math::var> polar<stan::math::var>(const stan::math::var& r,
                                                const stan::math::var& theta) {
  if (!(r >= 0) || stan::math::is_inf(theta)) {
    return {numeric_limits<double>::quiet_NaN()};
  }
  return {r * cos(theta), r * sin(theta)};
}

/**
 * Return the natural exponent of the specified complex number.
 *
 * @param[in] complex argument
 * @return exponential of argument
 */
template <>
complex<stan::math::var> exp<stan::math::var>(
    const complex<stan::math::var>& z) {
  using stan::math::is_inf;
  using stan::math::is_nan;
  if (is_inf(z.real()) && z.real() > 0) {
    if (is_nan(z.imag()) || z.imag() == 0) {
      // (+inf, nan), (+inf, 0)
      return z;
    } else if (is_inf(z.imag()) && z.imag() > 0) {
      // (+inf, +inf)
      return {z.real(), numeric_limits<double>::quiet_NaN()};
    } else if (is_inf(z.imag()) && z.imag() < 0) {
      // (+inf, -inf)
      return {numeric_limits<double>::quiet_NaN(),
              numeric_limits<double>::quiet_NaN()};
    }
  }
  if (is_inf(z.real()) && z.real() < 0
      && (is_nan(z.imag()) || is_inf(z.imag()))) {
    // (-inf, nan), (-inf, -inf), (-inf, inf)
    return {0, 0};
  }
  if (is_nan(z.real()) && z.imag() == -0.0) {
    // (nan, -0)
    return z;
  }
  stan::math::var exp_re = exp(z.real());
  return {exp_re * cos(z.imag()), exp_re * sin(z.imag())};
}

/**
 * Return the natural logarithm of the specified complex value with a
 * branch cut along the negative real axis.
 *
 * @param z complex argument
 * @return natural logarithm of argument
 */
template <>
complex<stan::math::var> log<stan::math::var>(
    const complex<stan::math::var>& z) {
  using stan::math::is_inf;
  using stan::math::is_nan;
  double inf = numeric_limits<double>::infinity();
  double nan = numeric_limits<double>::quiet_NaN();
  if ((is_nan(z.real()) && is_inf(z.imag()))
      || (is_inf(z.real()) && is_nan(z.imag()))) {
    return {inf, nan};
  }
  stan::math::var r = sqrt(norm(z));
  stan::math::var theta = arg(z);
  return {log(r), theta};
}

/**
 * Return the base 10 logarithm of the specified complex value with a
 * branch cut along the negative real axis.
 *
 * @param z complex argument
 * @return base 10 logarithm of argument
 */
template <>
complex<stan::math::var> log10<stan::math::var>(
    const complex<stan::math::var>& z) {
  static const double inv_log_10 = 1 / log(10);
  // TODO(carpenter): remove upcast to var after * overloaded for double
  return log(z) * stan::math::var(inv_log_10);
}

template <>
complex<stan::math::var> pow<stan::math::var>(
    const complex<stan::math::var>& x, const complex<stan::math::var>& y) {
  return exp(y * log(x));
}

// Following std:: templates cannot be specialized in both g++ and clang++

// can be specialized in g++, but not in clang++
// complex<var> pow(const complex<var>& x, const var& y);
// complex<var> pow(const var& x, const complex<var>& y);
// complex<var> pow(const complex<var>& x, int y);

// can't be specialized in g++ or clang++
// complex<var> pow(const complex<var>& x, const complex<double>& y);
// complex<var> pow(const complex<double>& x, const complex<var>& y);
// complex<var> pow(const complex<var>& x, const double& y);
// complex<var> pow(const double& x, const complex<var>& y);

/**
 * Return the square root of the specified complex number with a
 * branch cut along the negative real axis.
 *
 * @param z complex number argument
 * @return square root of argument
 */
template <>
complex<stan::math::var> sqrt<stan::math::var>(
    const complex<stan::math::var>& z) {
  auto m = sqrt(hypot(z.real(), z.imag()));
  auto at = 0.5 * atan2(z.imag(), z.real());
  return {m * cos(at), m * sin(at)};
}

/**
 * Return the complex hyperbolic sine of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic sine of argument
 */
template <>
complex<stan::math::var> sinh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return stan::math::var(0.5) * (exp(z) - exp(-z));
}

/**
 * Return the complex hyperbolic cosine of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic cosine of argument
 */
template <>
complex<stan::math::var> cosh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return stan::math::var(0.5) * (exp(z) + exp(-z));
}

/**
 * Return the complex hyperbolic tangent of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic tangent of argument
 */
template <>
complex<stan::math::var> tanh<stan::math::var>(
    const complex<stan::math::var>& z) {
  auto exp_z = exp(z);
  auto exp_neg_z = exp(-z);
  return (exp_z - exp_neg_z) / (exp_z + exp_neg_z);
}

/**
 * Return the complex arc hyperbolic sine of the specified complex
 * argument with branch cuts outside the interval `[-i, i]` along the
 * imaginary axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic sine of argument
 */
template <>
complex<stan::math::var> asinh<stan::math::var>(
    const complex<stan::math::var>& z) {
  complex<double> y_d = asinh(stan::math::value_of(z));
  auto y = log(z + sqrt(stan::math::var(1) + z * z));
  return copysign(y, y_d);
}

/**
 * Return the complex arc hyperbolic cosine of the specified complex
 * argument with branch cuts at values less than 1 along the real
 * axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic cosine of argument
 */
template <>
complex<stan::math::var> acosh<stan::math::var>(
    const complex<stan::math::var>& z) {
  complex<double> y_d = acosh(stan::math::value_of(z));
  auto y = log(z + sqrt(z * z - stan::math::var(1)));
  return copysign(y, y_d);
}

/**
 * Return the complex arc hyperbolic tangent of the specified complex
 * argument with branch cuts outside the interval `[-i, i]` along the
 * real axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic tanget of argument
 */
template <>
complex<stan::math::var> atanh<stan::math::var>(
    const complex<stan::math::var>& z) {
  complex<double> y_d = atanh(stan::math::value_of(z));
  stan::math::var one{1};
  auto y = stan::math::var(0.5) * (log(one + z) - log(one - z));
  return copysign(y, y_d);
}

/**
 * Return the complex sine of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex sine of argument
 */
template <>
complex<stan::math::var> sin<stan::math::var>(
    const complex<stan::math::var>& z) {
  return neg_i_times(sinh(i_times(z)));
}

/**
 * Return the complex cosine of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex cosine of argument
 */
template <>
complex<stan::math::var> cos<stan::math::var>(
    const complex<stan::math::var>& z) {
  return cosh(i_times(z));
}

/**
 * Return the complex tangent of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex tangent of argument
 */
template <>
complex<stan::math::var> tan<stan::math::var>(
    const complex<stan::math::var>& z) {
  return neg_i_times(tanh(i_times(z)));
}

/**
 * Return the complex arc sine of the specified complex argument, with
 * branch cuts outside the interval `[-1, 1]` along the real axis.
 *
 * @param[in] z complex argument
 * @return complex arc sine of the argument
 */
template <>
complex<stan::math::var> asin<stan::math::var>(
    const complex<stan::math::var>& z) {
  complex<double> y_d = asin(stan::math::value_of(z));
  auto y = neg_i_times(asinh(i_times(z)));
  return copysign(y, y_d);
}

/**
 * Return the complex arc cosine of the specified complex argument, with
 * branch cuts outside the interval `[-1, 1]` along the real axis.
 *
 * @param[in] z complex argument
 * @return complex arc cosine of the argument
 */
template <>
complex<stan::math::var> acos<stan::math::var>(
    const complex<stan::math::var>& z) {
  return stan::math::var(0.5 * stan::math::pi()) - asin(z);
}

/**
 * Return the complex arc tangent of the specified complex argument,
 * with branch cuts outside the interval `[-i, i]` along the imaginary
 * axis.
 *
 * @param[in] z complex argument
 * @return complex arc tangent of the argument
 */
template <>
complex<stan::math::var> atan<stan::math::var>(
    const complex<stan::math::var>& z) {
  return neg_i_times(atanh(i_times(z)));
}

}  // namespace std

// SPECIALIZATIONS FOR ADL IN STAN

namespace stan {
namespace math {
// no partial fun specializations allowed in C++1y
// instead, overload and rely in argument-dependent lookup in stan::math
/**
 * Writes specified complex number to specified ostream in
 * form `(real, imag)`.
 *
 * @param o[in,out] character output stream
 * @param x[in] complex number to be inserted
 * @return specified output stream
 */
template <class CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<<(
    std::basic_ostream<CharT, Traits>& o,
    const std::complex<stan::math::var>& x) {
  std::basic_ostringstream<CharT, Traits> s;
  s.flags(o.flags());
  s.imbue(o.getloc());
  s.precision(o.precision());
  s << '(' << x.real() << "," << x.imag() << ')';
  return o << s.str();
}

template <class CharT, class Traits>
std::basic_istream<CharT, Traits>& operator>>(
    std::basic_istream<CharT, Traits>& is, std::complex<stan::math::var>& x) {
  std::complex<double> y;
  is >> y;
  x = y;
  return is;
}
}  // namespace math
}  // namespace stan

namespace stan {
namespace math {

/**
 * Return the sum of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return sum of the arguments
 */
std::complex<var> operator+(const var& lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const var& lhs, const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const std::complex<var>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const std::complex<var>& lhs,
                            const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const std::complex<var>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const std::complex<double>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
std::complex<var> operator+(const std::complex<double>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator+(const std::complex<var>& lhs, T rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator+(T lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}

/**
 * Return the difference of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return difference of the arguments
 */
std::complex<var> operator-(const var& lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const var& lhs, const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const std::complex<var>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const std::complex<var>& lhs,
                            const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const std::complex<var>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const std::complex<double>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
std::complex<var> operator-(const std::complex<double>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator-(const std::complex<var>& lhs, T rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator-(T lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y -= rhs;
  return y;
}

/**
 * Return the product of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return product of the arguments
 */
std::complex<var> operator*(const var& lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const var& lhs, const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const std::complex<var>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const std::complex<var>& lhs,
                            const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const std::complex<var>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const std::complex<double>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
std::complex<var> operator*(const std::complex<double>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator*(const std::complex<var>& lhs, T rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator*(T lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y *= rhs;
  return y;
}

/**
 * Return the quotient of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return quotient of the arguments
 */
std::complex<var> operator/(const var& lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const var& lhs, const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const std::complex<var>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const std::complex<var>& lhs,
                            const std::complex<double>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const std::complex<var>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const std::complex<double>& lhs,
                            const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
std::complex<var> operator/(const std::complex<double>& lhs, const var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator/(const std::complex<var>& lhs, T rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> operator/(T lhs, const std::complex<var>& rhs) {
  std::complex<stan::math::var> y(lhs);
  y /= rhs;
  return y;
}

/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument raised to the power of the second
 */
std::complex<var> pow(const std::complex<var>& x, const var& y) {
  return exp(y * log(x));
}
std::complex<var> pow(const var& x, const std::complex<var>& y) {
  return exp(y * log(x));
}
std::complex<var> pow(const std::complex<var>& x, int y) {
  return exp(y * log(x));
}

std::complex<var> pow(const std::complex<var>& x,
                      const std::complex<double>& y) {
  return exp(y * log(x));
}
std::complex<var> pow(const std::complex<double>& x,
                      const std::complex<var>& y) {
  return exp(y * log(x));
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> pow(const std::complex<var>& x, const T& y) {
  return exp(static_cast<double>(y) * log(x));
}
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
std::complex<var> pow(const T& x, const std::complex<var>& y) {
  return exp(y * log(static_cast<double>(x)));
}

}  // namespace math
}  // namespace stan
#endif
