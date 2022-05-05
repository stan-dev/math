#ifndef STAN_MATH_REV_FUN_LOG1P_HPP
#define STAN_MATH_REV_FUN_LOG1P_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1p.hpp>

namespace stan {
namespace math {

/**
 * The log (1 + x) function for variables (C99).
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \log (1 + x) = \frac{1}{1 + x}\f$.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 * @return The log of 1 plus the variable.
 */
inline auto log1p(const var a) {
  return make_callback_var(log1p(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (1.0 + a.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log1p(const T& a) {
  auto a_arena = to_arena(a);
  return make_callback_rev_matrix<T>(log1p(a_arena.val()), [a_arena](auto&& vi) mutable {
    a_arena.adj().array() += vi.adj().array() / (1.0 + a_arena.val().array());
  });
}

}  // namespace math
}  // namespace stan
#endif
