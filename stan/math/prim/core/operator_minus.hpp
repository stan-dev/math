#ifndef STAN_MATH_PRIM_CORE_OPERATOR_MINUS_HPP
#define STAN_MATH_PRIM_CORE_OPERATOR_MINUS_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

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
inline std::complex<U> operator-(const std::complex<U>& x) {
  return {-x.real(), -x.imag()};
}

}  // namespace math
}  // namespace stan

#endif
