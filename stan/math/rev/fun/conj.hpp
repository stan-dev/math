#ifndef STAN_MATH_REV_FUN_CONJ_HPP
#define STAN_MATH_REV_FUN_CONJ_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/conj.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the complex conjugate of the complex argument.
 *
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
inline std::complex<var> conj(const std::complex<var>& z) {
  return internal::complex_conj(z);
}

}  // namespace math
}  // namespace stan

#endif
