#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

// doesn't need namespace std

#include <stan/math/rev/core.hpp>
#include <cmath>
#include <complex>
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
   * Adds other to real part of this.
   *
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  complex<stan::math::var>& operator+=(const stan::math::var& other) {
    re_ += other;
    return *this;
  }

  complex<stan::math::var>& operator-=(const stan::math::var& other) {
    re_ -= other;
    return *this;
  }

  complex<stan::math::var>& operator*=(const stan::math::var& other) {
    re_ *= other;
    return *this;
  }

  complex<stan::math::var>& operator/=(const stan::math::var& other) {
    re_ /= other;
    return *this;
  }

  // // (5)
  // template <typename T>
  // complex<stan::math::var>& operator+=(const std::complex<T>& other) {
  //   re_ += other.re_;
  //   im_ += other.im_;
  //   return *this;
  // }
  // // (6)
  // template <typename T>
  // complex<stan::math::var>& operator-=(const std::complex<T>& other) {
  //   re_ -= other.re_;
  //   im_ -= other.im_;
  //   return *this;
  // }
  // // (7)
  // template <typename T>
  // complex<stan::math::var>& operator*=(const std::complex<T>& other) {
  //   re_ = re_ * other.re_ - im_ * other.im_;
  //   im_ = re_ * other.im_ + other.re_ * im_;
  //   return *this;
  // }
  // // (8)
  // template <typename T>
  // complex<stan::math::var>& operator/=(const std::complex<T>& other) {
  //   stan::math::var sum_sq_im = square(im_) + square(other.im_);
  //   re_ = (re_ * other.re_ + im_ * other.im_) / sum_sq_im;
  //   im_ = (im_ * other.re_ - re_ * other.im_) / sum_sq_im;
  //   return *this;
  // }
};

// // non-member function specializations

// // (1)
// template <>
// std::complex<stan::math::var>
// operator+(const std::complex<stan::math::var>& val) {
//   return val;
// }
// // (2)
// template <>
// std::complex<stan::math::var>
// operator-(const std::complex<stan::math::var>& val) {
//   return std::complex<stan::math::var>(-val.re_, -val.im_);
// }

// // (1)
// template <>
// std::complex<stan::math::var>
// operator+(const std::complex<stan::math::var>& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y += rhs;
//   return y;
// }
// // (2)
// template <>
// std::complex<stan::math::var>
// operator+(const std::complex<stan::math::var>& lhs,
//           const stan::math::var& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y += rhs;
//   return y;
// }
// // (3)
// template <>
// std::complex<stan::math::var>
// operator+(const stan::math::var& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y += rhs;
//   return y;
// }
// // (4)
// template <>
// std::complex<stan::math::var>
// operator-(const std::complex<stan::math::var>& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y -= rhs;
//   return y;
// }
// // (5)
// template <>
// std::complex<stan::math::var>
// operator-(const std::complex<stan::math::var>& lhs,
//           const stan::math::var& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y -= rhs;
//   return y;
// }
// // (6)
// template <>
// std::complex<stan::math::var>
// operator-(const stan::math::var& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y -= rhs;
//   return y;
// }
// // (7)
// template <>
// std::complex<stan::math::var>
// operator*(const std::complex<stan::math::var>& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y *= rhs;
//   return y;
// }
// // (8)
// template <>
// std::complex<stan::math::var>
// operator*(const std::complex<stan::math::var>& lhs,
//           const stan::math::var& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y *= rhs;
//   return y;
// }
// // (9)
// template <>
// std::complex<stan::math::var>
// operator/(const stan::math::var& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y /= rhs;
//   return y;
// }
// // (10)
// template <>
// std::complex<stan::math::var>
// operator/(const std::complex<stan::math::var>& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y /= rhs;
//   return y;
// }
// // (11)
// template <>
// std::complex<stan::math::var>
// operator/(const std::complex<stan::math::var>& lhs,
//           const stan::math::var& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y /= rhs;
//   return y;
// }
// // (12)
// template <>
// std::complex<stan::math::var>
// operator/(const stan::math::var& lhs,
//           const std::complex<stan::math::var>& rhs) {
//   std::complex<stan::math::var> y(lhs);
//   y /= rhs;
//   return y;
// }

// // (1)
// template <>
// bool operator==(const std::complex<stan::math::var>& lhs,
//                 const std::complex<stan::math::var>& rhs) {
//   return lhs.re_ == rhs.re_ && lhs.im_ == rhs.im_;
// }
// // (2)
// template <>
// bool operator==(const std::complex<stan::math::var>& lhs,
//                 const stan::math::var& rhs) {
//   return lhs.re_ == rhs && lhs.im_ == 0;
// }
// // (3)
// template <>
// bool operator==(const stan::math::var& lhs,
//                 const std::complex<stan::math::var>& rhs) {
//   return !(lhs == rhs.re_ && 0 == rhs.im_);
// }
// // (4)
// template <>
// bool operator!=(const std::complex<stan::math::var>& lhs,
//                 const std::complex<stan::math::var>& rhs) {
//   return !(lhs.re_ == rhs.re_ && lhs.im_ == rhs.im_);
// }
// // (5)
// template <>
// bool operator!=(const std::complex<stan::math::var>& lhs,
//                 const stan::math::var& rhs) {
//   return !(lhs.re_ == rhs && lhs.im_ == 0);
// }
// // (6)
// template <>
// bool operator!=(const stan::math::var& lhs,
//                 const std::complex<stan::math::var>& rhs) {
//   return !(lhs == rhs.re_ && 0 == rhs.im_);
// }

// // (1)
// template<class CharT, class Traits>
// basic_ostream<CharT, Traits>&
// operator<< <stan::math::var, CharT, Traits>(basic_ostream<charT, traits>& o,
//                                             const complex<T>& x) {
//     basic_ostringstream<CharT, Traits> s;
//     s.flags(o.flags());
//     s.imbue(o.getloc());
//     s.precision(o.precision());
//     s << '(' << x.real() << "," << x.imag() << ')';
//     return o << s.str();
// }
// // (2)
// template <class CharT, class Traits>
// std::basic_istream<CharT, Traits>&
// operator>> <stan::math::var, CharT, Traits>(
//     std::basic_istream<CharT, Traits>& is,
//     std::complex<stan::math::var>& x) {
//   // TODO(carpenter): real implementation here
//   return is;
// }

// // (1)
// template <>
// stan::math::var real<stan::math::var>(const std::complex<stan::math::var>& z)
// {
//   return z.re_;
// }

// // (1)
// template <>
// stan::math::var imag<stan::math::var>(const std::complex<stan::math::var>& z)
// {
//   return z.im_;
// }

// // (1)
// template <>
// stan::math::var abs<stan::math::var>(const std::complex<stan::math::var>& z)
// {
//   using std::hypot;
//   return hypot(std::real(z), std::imag(z));
// }

// // (1)
// template <>
// stan::math::var arg<stan::math::var>(const std::complex<stan::math::var>& z)
// {
//   using std::atan2;
//   return atan2(std::imag(z), std::real(z));
// }

// // (1)
// template <>
// stan::math::var norm<stan::math::var>(const std::complex<stan::math::var>& z)
// {
//   return square(std::real(z)) + square(std::imag(z));
// }

// // (1)
// template <>
// std::complex<stan::math::var> conj<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   return {z.re_, -z.im_};
// }

// // (1)
// template <>
// std::complex<stan::math::var> proj<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   using std::isinf;
//   using std::copysign;
//   if (isinf(z.re_) || isinf(z.im_)) {
//     return { std::numeric_limits<stan::math::var>::infinity(),
//           copysign(0.0, z.im_) };  // (magnitude, sign)
//   }
//   return z;
// }

// // TODO(carpenter): add mixed forms with primitives for efficiency
// // (1)
// template <>
// std::complex<stan::math::var> polar<stan::math::var>(const stan::math::var&
// r,
//                                              const stan::math::var& theta) {
//   using std::cos;
//   using std::sin;
//   if (!(r >= 0)) {
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//   return {r * cos(theta), r * sin(theta)};
// }

// // (1)
// template <>
// std::complex<stan::math::var> exp<stan::math::var>(
//     const std::complex<stan::math::var>& z) {
//   using std::exp;
//   stan::math::var exp_re = exp(z.re_);
//   return {exp_re * cos(z.im_), exp_re * sin(z.im_)};
// }

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
//   stan::math::var hypt = hypot(z.re_, z.im_);
//   stan::math::var sqrt_hypot = sqrt(hypt);
//   stan::math::var half_atan2_z = 0.5 * atan2(z.re_, z.im_);
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
