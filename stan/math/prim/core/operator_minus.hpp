#ifndef STAN_MATH_PRIM_CORE_OPERATOR_MINUS_HPP
#define STAN_MATH_PRIM_CORE_OPERATOR_MINUS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

/**
 * Return the negation of the argument.
 *
 * @tparam U value type argument
 * @param x argument
 * @return negation of the argument
 */
template <typename U, require_autodiff_t<U>>
inline stan::math::complex<U> operator-(const stan::math::complex<U>& x) {
  return {-x.real(), -x.imag()};
}

}  // namespace math
}  // namespace stan

#endif
