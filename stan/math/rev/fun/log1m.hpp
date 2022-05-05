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
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 * @return The variable representing log of 1 minus the variable.
 */
inline auto log1m(const var a) {
  return make_callback_var(log1m(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (a.val() - 1.0);
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log1m(const T& a) {
  auto a_arena = to_arena(a);
  return make_callback_rev_matrix<T>(
      log1m(a_arena.val()), [a_arena](auto&& vi) mutable {
        a_arena.adj().array()
            += vi.adj().array() / (a_arena.val().array() - 1.0);
      });
}

}  // namespace math
}  // namespace stan
#endif
