#ifndef STAN_MATH_REV_CPLX_OPERATOR_DIVISION_HPP
#define STAN_MATH_REV_CPLX_OPERATOR_DIVISION_HPP

#include <stan/math/rev/cplx/var.hpp>
#include <stan/math/rev/scal/fun/pow.hpp>
#include <stan/math/prim/cplx/operator_division.hpp>

namespace stan {
namespace math {

/**
 * Return the result of dividing the first argument by the second.
 * This overload exists because libc++ uses logb and scalbn for
 * complex division, which aren't defined for fvars.
 *
 * @tparam T type of complex var
 * @param t first argument
 * @param u second argument
 * @return first argument divided by second argument
 */
inline std::complex<var> operator/(
    std::complex<var> const& t,
    std::complex<var> const& u) {
  return stan::math::operator_division(t, u);  // no recursion
}

}  // namespace math
}  // namespace stan
#endif
