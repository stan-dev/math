#ifndef STAN_MATH_FWD_FUN_ARG_HPP
#define STAN_MATH_FWD_FUN_ARG_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/arg.hpp>
#include <complex>

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
inline fvar<T> arg(const std::complex<fvar<T>>& z) {
  return internal::complex_arg(z);
}

}  // namespace math
}  // namespace stan

#endif
