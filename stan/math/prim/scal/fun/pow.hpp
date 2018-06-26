#ifndef STAN_MATH_PRIM_SCAL_FUN_POW_HPP
#define STAN_MATH_PRIM_SCAL_FUN_POW_HPP

#include <stan/math/prim/scal/fun/complex_promote.hpp>
#include <stan/math/prim/scal/meta/is_complex.hpp>
#include <stan/math/prim/scal/meta/rm_complex.hpp>
#include <cmath>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

template <class T, class U,
          std::enable_if_t<!std::is_same<rm_complex_t<T>,
                                         rm_complex_t<U>>::value>* = nullptr>
inline auto pow(std::complex<T> const& t, U const& u) {
  return pow(complex_promote<U>(t), complex_promote<T>(u));
}

template <
    class T, class U,
    std::enable_if_t<!std::is_same<rm_complex_t<T>, rm_complex_t<U>>::value
                     && !is_complex<T>::value>* = nullptr>
inline auto pow(T const& t, std::complex<U> const& u) {
  return pow(complex_promote<U>(t), complex_promote<T>(u));
}

}  // namespace math
}  // namespace stan
#endif
