#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/std_numeric_limits.hpp>
#include <stan/math/rev/core/std_isnan.hpp>
#include <stan/math/rev/core/std_isinf.hpp>
#include <complex>

namespace std {

/**
 * Specialization of the std::complex<stan::math::var>
 * constructor.
 *
 * The base template implementation will create uninitialized
 * stan::math::vars for the default arguments by calling the empty
 * constructor. This needs to be specialized because uninitialized
 * vars can cause seg faults.
 *
 * Promotion from primitives to stan::math::var applies and this
 * constructor will handle the primitivies.
 *
 * @param re real component
 * @param im imaginary component
 */
template <>
constexpr complex<stan::math::var>::complex(const stan::math::var& re,
                                            const stan::math::var& im) {
  real(stan::math::is_uninitialized(re) ? stan::math::var{0.0} : re);
  imag(stan::math::is_uninitialized(im) ? stan::math::var{0.0} : im);
}

/**
 * Specialization of operator= for stan::math::var.
 *
 * The base template implementation will leave the imaginary component
 * as an uninitialized stan::math::var, which can cause seg faults.
 * This implementation fixes that problem by forcing the imaginary
 * component to be initialized to 0.0.
 *
 *
 */
template <>
complex<stan::math::var>& complex<stan::math::var>::operator=(
    const stan::math::var& c) {
  real(c);
  imag(0.0);
  return *this;
}

template <>
inline stan::math::var norm(const complex<stan::math::var>& c) {
  return stan::math::square(c.real()) + stan::math::square(c.imag());
}

template <>
complex<stan::math::var> operator/(const complex<stan::math::var>& z,
                                   const complex<stan::math::var>& w) {
  return complex<stan::math::var>{z.real() * w.real() + z.imag() * w.imag(),
                                  z.imag() * w.real() - z.real() * w.imag()}
         / norm(w);
}

template <>
complex<stan::math::var> operator*(const complex<stan::math::var>& z,
                                   const complex<stan::math::var>& w) {
  return complex<stan::math::var>{z.real() * w.real() - z.imag() * w.imag(),
                                  z.real() * w.imag() + z.imag() * w.real()};
}

template <>
template <>
complex<stan::math::var>& complex<stan::math::var>::operator/=(
    const complex<stan::math::var>& y) {
  *this = *this / y;
  return *this;
}

template <>
template <>
complex<stan::math::var>& complex<stan::math::var>::operator*=(
    const complex<stan::math::var>& y) {
  *this = *this * y;
  return *this;
}

template <>
inline bool operator==(const complex<stan::math::var>& x,
                       const complex<stan::math::var>& y) {
  return x.real() == y.real() && x.imag() == y.imag();
  ;
}

template <>
inline bool operator==(const complex<stan::math::var>& x,
                       const stan::math::var& y) {
  return x.real() == y && x.imag() == 0;
}

template <>
inline bool operator==(const stan::math::var& x,
                       const complex<stan::math::var>& y) {
  return y == x;
}

inline bool operator==(const complex<stan::math::var>& x, double y) {
  return x.real() == y && x.imag() == 0;
}

inline bool operator==(double x, const complex<stan::math::var>& y) {
  return y == x;
}

template <>
inline bool operator!=(const complex<stan::math::var>& x,
                       const complex<stan::math::var>& y) {
  return x.real() != y.real() || x.imag() != y.imag();
}

inline bool operator!=(const complex<stan::math::var>& x,
                       const stan::math::var& y) {
  return (x.real() != y || x.imag() != 0);
}

inline bool operator!=(const stan::math::var& x,
                       const complex<stan::math::var>& y) {
  return y != x;
}

inline bool operator!=(const complex<stan::math::var>& x, double y) {
  return (x.real() != y || x.imag() != 0);
  ;
}

inline bool operator!=(double x, const complex<stan::math::var>& y) {
  return y != x;
}

inline int isinf(const std::complex<stan::math::var>& a) {
  return stan::math::is_inf(a.real().val());
}

inline int isnan(const std::complex<stan::math::var>& a) {
  return stan::math::is_nan(a.real().val());
}

template <>
inline stan::math::var abs(const complex<stan::math::var>& c) {
  return stan::math::sqrt(norm(c));
}

template <>
inline complex<stan::math::var> conj(const complex<stan::math::var>& c) {
  return complex<stan::math::var>(c.real(), -c.imag());
}

template <>
inline complex<stan::math::var> proj(const complex<stan::math::var>& c) {
  if (isinf(c.real()) || isinf(c.imag()))
    return std::complex<stan::math::var>(stan::math::positive_infinity(),
                                         copysign(0.0, c.imag().val()));
  return c;
}

}  // namespace std
#endif
