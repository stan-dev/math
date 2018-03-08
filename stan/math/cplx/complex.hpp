#ifndef STAN_MATH_CPLX_COMPLEX_HPP
#define STAN_MATH_CPLX_COMPLEX_HPP

#include <complex>
#include <type_traits>

// this file applies to both forward and reverse mode, but isn't under mix
// because it should be included from var.hpp and fvar.hpp

namespace stan {
namespace math {
namespace internal {

template <class T>
class complex;  // forward declaration of stan's complex

template <class>
struct is_cplx : std::false_type {};
template <class T>
struct is_cplx<std::complex<T>> : std::true_type {};
template <class T>
struct is_cplx<complex<T>> : std::true_type {};

/// This class exists purely to forward the interface of std::complex and serve
/// as a tag for the free functions below. Without this class, the free
/// functions  would have to be written to extend std::complex, which would
/// change its semantics, or would have to be specialized to only stan types.
template <class T>
struct complex : std::complex<T> {
  using std::complex<T>::complex;  ///< inherit all std::complex ctors
  complex(const std::complex<T>& t)
      : std::complex<T>(t){};  ///< allow promotion
  /**without this converting ctor, copy initialization of std::complex<z_(f)var>
   * from an (f)var fails, because the default implementation in the
   * std::complex template requires the argument type to match the template
   * type. This is the same reason why std::complex<double> can't be multiplied
   * by an int.
   * tparam R type of real argument. Won't fire if complex.
   * tparam I type of imaginary argument.
   * param[in] real, the pure real component of the complex number.
   * param[in] imag, the pure imaginary component of the complex number*/
  template <class R = T, class I = T,
            std::enable_if_t<!is_cplx<R>::value>* = nullptr>
  complex(const R real = 0, const I imag = 0) : std::complex<T>(real, imag) {}
};

// The free functions that follow extend stan::math::internal::complex beyond
// the  interface of std::complex. This is necessary because the std::complex
// template from which stan inherits does not allow simple operations such as
// std::complex<double>()*2, because 2 is an int, and the functions that handle
// std::complex only work if the  template type T (double in this case)
// appears in both function arguments. The compiler doesn't promote the int to
// double because the template type isn't found in time for that to happen.

// Since stan's complex is working with vars and similar classes, it would also
// suffer the same problem when the other argument is anything except a var
// without extending the binary functions as we have below. In most cases,
// the functions just delegate to single argument versions in std::complex.

template <class T, class U>
inline complex<T> operator+(complex<T> t, U u) {
  return t += u;
}

template <class T, class U>
inline complex<T> operator+(U u, complex<T> t) {
  return t += u;
}

template <class T, class U>
inline complex<T> operator*(complex<T> t, U u) {
  return t *= u;
}

template <class T, class U>
inline complex<T> operator*(U u, complex<T> t) {
  return t *= u;
}

template <class T, class U>
inline complex<T> operator-(complex<T> t, U u) {
  return t -= u;
}

template <class T, class U>
inline complex<T> operator-(U u, complex<T> t) {
  return (t -= u) * (-1);
}

template <class T, class U>
inline complex<T> operator/(complex<T> t, U u) {
  return t /= u;
}

template <class T, class U>
inline complex<T> operator/(U u, complex<T> t) {
  return complex<T>(u) /= t;
}

template <class T, class U>
inline bool operator==(complex<T> t, complex<U> u) {
  return t.real() == u.real() && t.imag() == u.imag();
}

template <class T, class U>
inline bool operator==(std::complex<U> u, complex<T> t) {
  return complex<T>(u) == t;
}

template <class T, class U>
inline bool operator==(complex<T> t, std::complex<U> u) {
  return complex<T>(u) == t;
}

template <class T, class U>
inline bool operator==(complex<T> t, U u) {
  return t.real() == u && t.imag() == T();
}

template <class T, class U>
inline bool operator==(U u, complex<T> t) {
  return t.real() == u && t.imag() == T();
}

template <class T, class U>
inline bool operator!=(complex<T> t, U u) {
  return !(t == u);
}

template <class T, class U>
inline bool operator!=(U u, complex<T> t) {
  return !(t == u);
}

template <class T, class U>
inline complex<T> pow(complex<T> t, complex<U> u) {
  return t == T() ? T() : std::exp(u * std::log(t));
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
