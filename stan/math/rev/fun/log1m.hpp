#ifndef STAN_MATH_REV_FUN_LOG1M_HPP
#define STAN_MATH_REV_FUN_LOG1M_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1m.hpp>

namespace stan {
namespace math {

/**
 * The log (1 - x) function for variables.
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \log (1 - x) = -\frac{1}{1 - x}\f$.
 *
 * @param a The variable.
 * @return The variable representing log of 1 minus the variable.
 */
inline var log1m(const var& a) {
  return make_callback_var(log1m(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (a.val() - 1);
  });
}

/**
 * Return the elementwise log of 1 - x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise log of 1 - x
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log1m(const T& x) {
  return make_callback_var(
      stan::math::log1m(x.val()), [x](const auto& vi) mutable {
        x.adj() += (vi.adj().array() / (x.val().array() - 1.0)).matrix();
      });
}

}  // namespace math
}  // namespace stan
#endif
