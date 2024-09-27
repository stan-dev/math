#ifndef STAN_MATH_FWD_FUN_CONJ_HPP
#define STAN_MATH_FWD_FUN_CONJ_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/conj.hpp>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

/**
 * Return the phase angle of the complex argument.
 *
 * @tparam T value type of autodiff variable
 * @param[in] z argument
 * @return phase angle of the argument
 */
template <typename T>
inline stan::math::complex<fvar<T>> conj(
    const stan::math::complex<fvar<T>>& z) {
  return internal::complex_conj(z);
}

}  // namespace math
}  // namespace stan

#endif
