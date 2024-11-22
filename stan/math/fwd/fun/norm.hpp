#ifndef STAN_MATH_FWD_FUN_NORM_HPP
#define STAN_MATH_FWD_FUN_NORM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/prim/fun/norm.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the squared magnitude of the complex argument.
 *
 * @tparam T value type of autodiff variable
 * @param[in] z argument
 * @return phase squared magnitude of the argument
 */
template <typename T>
inline fvar<T> norm(const std::complex<fvar<T>>& z) {
  return internal::complex_norm(z);
}

}  // namespace math
}  // namespace stan

#endif
