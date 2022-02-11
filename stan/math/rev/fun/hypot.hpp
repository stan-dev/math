#ifndef STAN_MATH_REV_FUN_HYPOT_HPP
#define STAN_MATH_REV_FUN_HYPOT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/hypot.hpp>

namespace stan {
namespace math {

/**
 * Returns the length of the hypotenuse of a right triangle
 * with sides of the specified lengths (C99).
 *
 * The partial derivatives are given by
 *
 * \f$\frac{\partial}{\partial x} \sqrt{x^2 + y^2} = \frac{x}{\sqrt{x^2 +
 * y^2}}\f$, and
 *
 * \f$\frac{\partial}{\partial y} \sqrt{x^2 + y^2} = \frac{y}{\sqrt{x^2 +
 * y^2}}\f$.
 *
 * @param[in] a Length of first side.
 * @param[in] b Length of second side.
 * @return Length of hypotenuse.
 */
template <typename T1, typename T2,
          require_any_var_t<T1, T2>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline var hypot(const T1& a, const T2& b) {
  auto args_tuple = std::make_tuple(a, b);
  auto val_fun = [&](auto&& x, auto&& y) { using std::hypot; return hypot(x, y); };
  auto grad_fun_tuple = std::make_tuple(
    [&](auto&& adj, auto&& x, auto&& y) { return adj * x / hypot(x, y); },
    [&](auto&& adj, auto&& x, auto&& y) { return adj * y / hypot(x, y); }
  );
  return user_gradients(args_tuple, val_fun, grad_fun_tuple);
}

}  // namespace math
}  // namespace stan
#endif
