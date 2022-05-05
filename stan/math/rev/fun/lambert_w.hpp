#ifndef STAN_MATH_REV_FUN_LAMBERT_W_HPP
#define STAN_MATH_REV_FUN_LAMBERT_W_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lambert_w.hpp>

namespace stan {
namespace math {

/**
 * Return the Lambert W function on W0 branch applied to the specified variable.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Variable argument.
 * @return the Lambert W function (W0 branch) applied to the specified argument.
 */
inline auto lambert_w0(const var a) {
  return make_callback_var(lambert_w0(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (a.val() + exp(vi.val()));
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto lambert_w0(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(
      lambert_w0(x_arena.val()), [x_arena](auto& vi) mutable {
        x_arena.adj().array()
            += vi.adj().array() / (x_arena.val() + exp(vi.val())).array();
      });
}

/**
 * Return the Lambert W function on W-1 branch applied to the specified
 * variable.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Variable argument.
 * @return the Lambert W function (W-1 branch) applied to the specified
 * argument.
 */
inline auto lambert_wm1(const var a) {
  return make_callback_var(lambert_wm1(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (a.val() + exp(vi.val()));
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto lambert_wm1(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(
      lambert_wm1(x_arena.val()), [x_arena](auto&& vi) mutable {
        x_arena.adj().array()
            += vi.adj().array() / (x_arena.val() + exp(vi.val())).array();
      });
}

}  // namespace math
}  // namespace stan

#endif
