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
 * @param a Variable argument.
 * @return the Lambert W function (W0 branch) applied to the specified argument.
 */
inline var lambert_w0(const var& a) {
  return make_callback_var(lambert_w0(a.val()), [a](auto& vi) mutable {
    a.adj() += (vi.adj() / (a.val() + exp(vi.val())));
  });
}

/**
 * Return the Lambert W function on W-1 branch applied to the specified
 * variable.
 *
 * @param a Variable argument.
 * @return the Lambert W function (W-1 branch) applied to the specified
 * argument.
 */
inline var lambert_wm1(const var& a) {
  return make_callback_var(lambert_wm1(a.val()), [a](auto& vi) mutable {
    a.adj() += (vi.adj() / (a.val() + exp(vi.val())));
  });
}

}  // namespace math
}  // namespace stan

#endif
