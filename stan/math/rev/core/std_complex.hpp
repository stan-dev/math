#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

// doesn't need namespace std

#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/atan2.hpp>
#include <stan/math/rev/scal/fun/cos.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/rev/scal/fun/sin.hpp>
#include <stan/math/rev/scal/fun/hypot.hpp>
#include <stan/math/rev/scal/fun/is_inf.hpp>
#include <stan/math/rev/scal/fun/is_nan.hpp>
#include <stan/math/rev/scal/fun/square.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>

#include <iostream>

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class std::complex<stan::math::var> {
  stan::math::var re_;
  stan::math::var im_;

 public:
  /**
   * Type of real and imaginary parts.
   */
  typedef stan::math::var value_type;

  /**
   * Constructs complex number from real and imaginary parts.
   *
   * @param[in] re the real part
   * @param[in] im the imaginary part
   */
  complex(const stan::math::var& re = stan::math::var(0),
          const stan::math::var& im = stan::math::var(0))
      : re_(re), im_(im) {}

  /**
   * Constructs complex number with the contents of other.
   *
   * @param[in] other another complex to use as source
   */
  complex(const std::complex<stan::math::var>& other)
      : re_(other.re_), im_(other.im_) {}

  /**
   * Constructs the complex number from the specified complex number
   * of a different type.
   *
   * @param[in] other another complex to use as source
   */
  template <typename T>
  complex(const std::complex<T>& other)
      : re_(other.real()), im_(other.imag()) {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  complex<stan::math::var>& operator=(const stan::math::var& x) {
    re_ = x;
    im_ = 0;
    return *this;
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @param[in] x complex value to assign
   * @return this complex number
   */
  complex<stan::math::var>& operator=(const complex<stan::math::var>& x) {
    re_ = x.re_;
    im_ = x.im_;
    return *this;
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam T value type of argument (must be assignable to
   * `stan::math::var`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename T>
  complex<stan::math::var>& operator=(const complex<T>& x) {
    re_ = x.real();
    im_ = x.imag();
    return *this;
  }

  /**
   * Return the real part.
   *
   * @return the real part
   */
  stan::math::var real() const { return re_; }

  /**
   * Set the real part to the specified value.
   *
   * @param[in] x the value to set the real part to
   */
  void real(const stan::math::var& x) { re_ = x; }

  /**
   * Return the imaginary part.
   *
   * @return the imaginary part
   */
  stan::math::var imag() const { return im_; }

  /**
   * Set the imaginary part to the specified value.
   *
   * @param[in] x the value to set the imaginary part to
   */
  void imag(const stan::math::var& x) { im_ = x; }

  /**
   * Adds other to this.
   *
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  complex<stan::math::var>& operator+=(const stan::math::var& other) {
    re_ += other;
    return *this;
  }

  /**
   * Subtracts other from this.
   *
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  complex<stan::math::var>& operator-=(const stan::math::var& other) {
    re_ -= other;
    return *this;
  }

  /**
   * Multiplies other by this.
   *
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  complex<stan::math::var>& operator*=(const stan::math::var& other) {
    re_ *= other;
    im_ *= other;
    return *this;
  }

  /**
   * Divides other by this.
   *
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  complex<stan::math::var>& operator/=(const stan::math::var& other) {
    re_ /= other;
    im_ /= other;
    return *this;
  }

  /**
   * Adds other to this.
   *
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename T>
  complex<stan::math::var>& operator+=(const std::complex<T>& other) {
    re_ += other.real();
    im_ += other.imag();
    return *this;
  }

  /**
   * Subtracts other from this.
   *
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename T>
  complex<stan::math::var>& operator-=(const std::complex<T>& other) {
    re_ -= other.real();
    im_ -= other.imag();
    return *this;
  }

  /**
   * Multiplies this by other.
   *
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename T>
  complex<stan::math::var>& operator*=(const std::complex<T>& other) {
    stan::math::var re_temp = re_ * other.real() - im_ * other.imag();
    im_ = re_ * other.imag() + other.real() * im_;
    re_ = re_temp;
    return *this;
  }

  /**
   * Divides this by other.
   *
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename T>
  complex<stan::math::var>& operator/=(const std::complex<T>& other) {
    using stan::math::square;
    stan::math::var sum_sq_im = square(other.real()) + square(other.imag());
    stan::math::var re_temp
        = (re_ * other.real() + im_ * other.imag()) / sum_sq_im;
    im_ = (im_ * other.real() - re_ * other.imag()) / sum_sq_im;
    re_ = re_temp;
    return *this;
  }
};

/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <>
struct std::iterator_traits<stan::math::var> {
  /**
   * Iterator category for traits.
   */
  typedef std::random_access_iterator_tag iterator_category;

  /**
   * Type for difference between pointers.
   */
  typedef std::ptrdiff_t difference_type;

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

/**
 * Return the value of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return a copy of the argument
 */
template <>
std::complex<stan::math::var> std::operator+<stan::math::var>(
    const std::complex<stan::math::var>& val) {
  return val;
}

/**
 * Return the negation of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return negated argument
 */
template <>
std::complex<stan::math::var> std::operator-<stan::math::var>(
    const std::complex<stan::math::var>& val) {
  return std::complex<stan::math::var>(-val.real(), -val.imag());
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
std::complex<stan::math::var> std::operator+<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator+<stan::math::var>(
    const std::complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  std::complex<stan::math::var> y(lhs);
  y += rhs;
  return y;
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs scalar first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
std::complex<stan::math::var> std::operator+<stan::math::var>(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator-<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator-<stan::math::var>(
    const std::complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator-<stan::math::var>(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator*<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator*<stan::math::var>(
    const std::complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator*<stan::math::var>(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator/<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator/<stan::math::var>(
    const std::complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  std::complex<stan::math::var> y(lhs);
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
std::complex<stan::math::var> std::operator/(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
  std::complex<stan::math::var> y(lhs);
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
bool std::operator==<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
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
bool std::operator==<stan::math::var>(const std::complex<stan::math::var>& lhs,
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
bool std::operator==<stan::math::var>(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
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
bool std::operator!=<stan::math::var>(
    const std::complex<stan::math::var>& lhs,
    const std::complex<stan::math::var>& rhs) {
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
bool std::operator!=<stan::math::var>(const std::complex<stan::math::var>& lhs,
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
bool std::operator!=<stan::math::var>(
    const stan::math::var& lhs, const std::complex<stan::math::var>& rhs) {
  return !(lhs == rhs.real() && 0 == rhs.imag());
}

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

/**
 * Return the real component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
stan::math::var std::real<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  return z.real();
}

/**
 * Return the imaginary component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
stan::math::var std::imag<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  return z.imag();
}

/**
 * Return the magnitude of the specified complex number.
 *
 * @param[in] complex argument
 * @return magnitude of argument
 */
template <>
stan::math::var std::abs<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  using std::hypot;
  return hypot(std::real(z), std::imag(z));
}

/**
 * Return the phase angle of the specified complex number.
 *
 * @param[in] z complex argument
 * @return phase angle of argument
 */
template <>
stan::math::var std::arg<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  using std::atan2;
  return atan2(std::imag(z), std::real(z));
}

/**
 * Return the squared magnitude of the specified complex number.
 *
 * @param[in] z complex argument
 * @return squared magnitude of argument
 */
template <>
stan::math::var std::norm<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  using stan::math::square;
  return square(std::real(z)) + square(std::imag(z));
}

/**
 * Return the complex conjugate of the specified complex number.
 *
 * @param[in] z complex argument
 * @return complex conjugate of argument
 */
template <>
std::complex<stan::math::var> std::conj<stan::math::var>(
    const std::complex<stan::math::var>& z) {
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
std::complex<stan::math::var> std::proj<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  if (stan::math::is_inf(z.real()) || stan::math::is_inf(z.imag())) {
    return {std::numeric_limits<stan::math::var>::infinity(),
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
std::complex<stan::math::var> std::polar<stan::math::var>(
    const stan::math::var& r, const stan::math::var& theta) {
  using std::cos;
  using std::sin;
  if (!(r >= 0) || stan::math::is_inf(theta)) {
    return {std::numeric_limits<double>::quiet_NaN()};
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
std::complex<stan::math::var> std::exp<stan::math::var>(
    const std::complex<stan::math::var>& z) {
  using stan::math::is_inf;
  using stan::math::is_nan;
  using std::exp;
  if (is_inf(z.real()) && z.real() > 0) {
    if (is_nan(z.imag()) || z.imag() == 0) {
      // (+inf, nan), (+inf, 0)
      return z;
    } else if (is_inf(z.imag())) {
      // (+inf, -inf),  (+inf, +inf)
      return {z.real(), std::numeric_limits<double>::quiet_NaN()};
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

// // (1)
// template <>
// std::complex<stan::math::var> log<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   using std::log;
//   using std::sqrt;
//   stan::math::var r = sqrt(norm(z));
//   stan::math::var theta = arg(z);
//   return {log(r), theta};
// }

// // (1)
// template <>
// std::complex<stan::math::var> log10<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   using std::log;
//   return log(z) / log(10);
// }

// // (1)
// template <>
// std::complex<stan::math::var> pow<stan::math::var>(
//     const std::complex<stan::math::var>& x,
//     const std::complex<stan::math::var>& y) {
//   using std::exp;
//   using std::log;
//   return exp(y * log(x));
// }
// // (2)
// template <>
// std::complex<stan::math::var> pow<stan::math::var>(
//     const std::complex<stan::math::var>& x, const stan::math::var& y) {
//   using std::exp;
//   using std::log;
//   return exp(y * log(x));
// }
// // (3)
// template <>
// std::complex<stan::math::var> pow<stan::math::var>(const stan::math::var& x,
//     const std::complex<stan::math::var>& y) {
//   using std::exp;
//   using std::log;
//   using stan::math::log;
//   return exp(y * log(x));
// }
// // (4)
// template <U>
// std::complex<stan::math::var> pow<stan::math::var, U>(
//     const std::complex<stan::math::var>& x,
//     const std::complex<U>& y) {
//   using std::exp;
//   using std::log;
//   return exp(y * log(x));
// }
// // (5)
// template <U>
// std::complex<stan::math::var> pow<stan::math::var, U>(
//     const std::complex<stan::math::var>& x, const U& y) {
//   using std::exp;
//   using std::log;
//   return exp(y * log(x));
// }
// // (6)
// template <T>
// std::complex<stan::math::var> pow<T, stan::math::var>(
//     const T& x, const std::complex<stan::math::var>& y) {
//   using std::exp;
//   using std::log;
//   using stan::math::log;
//   return exp(y * log(x));
// }

// // (1)
// template<>
// std::complex<stan::math::var> sqrt<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   using stan::math::sqrt;
//   using stan::math::hypot;
//   stan::math::var hypt = hypot(z.real(), z.imag());
//   stan::math::var sqrt_hypot = sqrt(hypt);
//   stan::math::var half_atan2_z = 0.5 * atan2(z.real(), z.imag());
//   return { sqrt_hypot * cos(half_atan2_z),
//         sqrt_hypot * sin(half_atan2_z) }
// }

// // (1)
// template <>
// std::complex<stan::math::var> sin<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return -i * std::sinh(i * z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> cos<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return std::cosh(i * z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> tan<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return -i * std::tanh(i * z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> asin<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return -i * std::asinh(i * z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> acos<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return 0.5 * stan::math::pi() - std::asin(z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> atan<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   static const std::complex<double> i{0, 1};
//   return -i * std::atanh(i * z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> sinh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return 0.5 * (std::exp(z) - std::exp(-z));
// }

// // (1)
// template <>
// std::complex<stan::math::var> cosh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return 0.5 * (std::exp(z) + std::exp(-z));
// }

// // (1)
// template <>
// std::complex<stan::math::var> tanh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   auto exp_z = std::exp(z);
//   auto exp_neg_z = std::exp(-z);
//   return (exp_z - exp_neg_z) / (exp_z + exp_neg_z);
// }

// // (1)
// template <>
// std::complex<stan::math::var> asinh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return std::log(z + std::sqrt(1 + z * z));
// }

// // (1)
// template <>
// std::complex<stan::math::var> acosh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return std::log(z + std::sqrt(z + 1) * std::sqrt(z - 1));
// }

// // (1)
// template <>
// std::complex<stan::math::var> atanh<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return 0.5 * (std::log(1 + z) - std::log(1 - z));
// }

#endif
