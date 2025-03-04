#ifndef STAN_MATH_REV_FUN_INV_PHI_HPP
#define STAN_MATH_REV_FUN_INV_PHI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv_Phi.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse of unit normal cumulative density function.
 *
 * The derivative is the reciprocal of unit normal density function,
 *
 * @param p Probability
 * @return The unit normal inverse cdf evaluated at p
 */
inline var inv_Phi(const var& p) {
  double val = inv_Phi(p.val());
  return make_callback_var(val, [p, val](auto& vi) mutable {
    p.adj() += vi.adj() * exp(-std_normal_lpdf(val));
  });
}

/**
 * Return the elementwise inverse of unit normal cumulative density function.
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param p Probability vector
 * @return Elementwise unit normal inverse cdf
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto inv_Phi(const T& p) {
  auto arena_rtn = to_arena(inv_Phi(p.val()));
  return make_callback_var(arena_rtn, [p, arena_rtn](auto& vi) mutable {
    auto deriv
        = arena_rtn.unaryExpr([](auto x) { return exp(-std_normal_lpdf(x)); });
    p.adj() += elt_multiply(vi.adj(), deriv);
  });
}

}  // namespace math
}  // namespace stan
#endif
