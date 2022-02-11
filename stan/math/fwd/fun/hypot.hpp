#ifndef STAN_MATH_FWD_FUN_HYPOT_HPP
#define STAN_MATH_FWD_FUN_HYPOT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/hypot.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the length of the hypotenuse of a right triangle with
 * opposite and adjacent side lengths given by the specified
 * arguments (C++11).  In symbols, if the arguments are
 * <code>1</code> and <code>x2</code>, the result is <code>sqrt(x1 *
 * x1 + x2 * x2)</code>.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return Length of hypotenuse of right triangle with opposite
 * and adjacent side lengths x1 and x2.
 */
template <typename T1, typename T2,
          require_any_fvar_t<T1, T2>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto hypot(const T1& a, const T2& b) {
  auto args_tuple = std::make_tuple(a, b);
  auto val_fun = [&](auto&& x, auto&& y) {using std::hypot; return hypot(x, y); };
  auto grad_fun_tuple = std::make_tuple(
    [&](auto&& x, auto&& y) { return x / hypot(x, y); },
    [&](auto&& x, auto&& y) { return y / hypot(x, y); }
  );
  return user_gradients(args_tuple, val_fun, grad_fun_tuple);
}

}  // namespace math
}  // namespace stan
#endif
