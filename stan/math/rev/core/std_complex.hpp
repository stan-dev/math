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
complex<stan::math::var> operator/(const complex<stan::math::var>& z,
                                   const complex<stan::math::var>& w) {
  return complex<stan::math::var>{z.real() / w.real(), z.imag() / w.imag()};
}

template <>
complex<stan::math::var> operator*(const complex<stan::math::var>& z,
                                   const complex<stan::math::var>& w) {
  return complex<stan::math::var>{z.real() * w.real(), z.imag() * w.imag()};
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

inline bool operator!=(const complex<stan::math::var>& x, double y) {
  return !(x == y);
}

inline bool operator!=(double x, const complex<stan::math::var>& y) {
  return !(x == y);
}

inline int isinf(const std::complex<stan::math::var>& a) {
  return stan::math::is_inf(a.real().val());
}

inline int isnan(const std::complex<stan::math::var>& a) {
  return stan::math::is_nan(a.real().val());
}

}  // namespace std
#endif
