#ifndef STAN_MATH_FWD_CPLX_OPERATOR_DIVISION_HPP
#define STAN_MATH_FWD_CPLX_OPERATOR_DIVISION_HPP

#include <stan/math/fwd/cplx/fvar.hpp>
#include <stan/math/fwd/scal/fun/pow.hpp>
#include <stan/math/prim/cplx/operator_division.hpp>

namespace stan {
namespace math {

/**
 * Return the result of dividing the first argument by the second.
 * This overload exists because libc++ uses logb and scalbn for
 * complex division, which aren't defined for fvars.
 *
 * @tparam T type of complex fvar value and tangent
 * @param t first argument
 * @param u second argument
 * @return first argument divided by second argument
 */
template <class T>
inline std::complex<fvar<T>> operator/(std::complex<fvar<T>> const& t,
                                       std::complex<fvar<T>> const& u) {
  return stan::math::operator_division(t, u);  // no recursion
}

}  // namespace math
}  // namespace stan
#endif
