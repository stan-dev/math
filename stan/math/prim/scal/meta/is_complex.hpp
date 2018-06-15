#ifndef STAN_MATH_PRIM_SCAL_META_IS_COMPLEX_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_COMPLEX_HPP

#include <complex>

namespace stan {
namespace math {

template <class T>
struct complex;  // forward declare stan's complex

}  // namespace math

/// trait to see if the template parameter is complex
template <class>
struct is_complex : 
 std::false_type {};

template <class T>
struct is_complex<std::complex<T>> 
 : std::true_type {};

template <class T>
struct is_complex<stan::math::complex<T>>
 : std::true_type {};

}  // namespace stan
#endif
