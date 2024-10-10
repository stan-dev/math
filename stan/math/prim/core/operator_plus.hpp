#ifndef STAN_MATH_PRIM_CORE_OPERATOR_PLUS_HPP
#define STAN_MATH_PRIM_CORE_OPERATOR_PLUS_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the argument.
 *
 * @tparam U value type argument
 * @param x argument
 * @return argument
 */
template <typename U, require_autodiff_t<U>>
inline stan::math::complex<U> operator+(const stan::math::complex<U>& x) {
  return x;
}

}  // namespace math
}  // namespace stan

#endif
