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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto lambert_w0(const var_value<T>& a) {
  return make_callback_var(lambert_w0(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj())
        += (as_array_or_scalar(vi.adj())
            / as_array_or_scalar(a.val() + exp(vi.val())));
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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto lambert_wm1(const var_value<T>& a) {
  return make_callback_var(lambert_wm1(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj())
        += (as_array_or_scalar(vi.adj())
            / as_array_or_scalar(a.val() + exp(vi.val())));
  });
}

}  // namespace math
}  // namespace stan

#endif
