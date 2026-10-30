#ifndef STAN_MATH_REV_FUNCTOR_SOLVE_ITP_HPP
#define STAN_MATH_REV_FUNCTOR_SOLVE_ITP_HPP

#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/solve_itp.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <tuple>

namespace stan {
namespace math {

/**
 * Reverse-mode overload of the ITP scalar root solver.
 *
 * The root is found with primitive values.  Gradients with respect to the
 * callback arguments are then supplied by the scalar implicit-function rule
 *
 *   dx / dtheta = -(partial f / partial theta) / (partial f / partial x)
 *
 * evaluated at the numerical root.  The bracket and all solver controls are
 * treated as data and are not differentiated.
 */
template <typename F, typename... Args, require_any_st_var<Args...>* = nullptr>
inline var solve_itp_tol(const F& f, double lower, double upper,
                         double x_tolerance, double kappa1, double kappa2,
                         int64_t n0, int64_t max_num_steps, std::ostream* msgs,
                         const Args&... args) {
  static constexpr const char* function = "solve_itp";

  auto arena_args_tuple = make_chainable_ptr(std::make_tuple(eval(args)...));
  auto args_vals_tuple = math::apply(
      [&](const auto&... arena_args) {
        return std::make_tuple(to_ref(value_of(arena_args))...);
      },
      *arena_args_tuple);

  const double root = math::apply(
      [&](const auto&... args_vals) {
        return solve_itp_tol(f, lower, upper, x_tolerance, kappa1, kappa2, n0,
                             max_num_steps, msgs, args_vals...);
      },
      args_vals_tuple);

  // Compute partial f / partial x with parameters held at their values.
  double f_x;
  {
    nested_rev_autodiff nested;
    var x_var = root;
    var residual = math::apply(
        [&](const auto&... args_vals) {
          return eval(f(x_var, msgs, args_vals...));
        },
        args_vals_tuple);
    residual.grad();
    f_x = x_var.adj();
  }

  check_finite(function, "partial derivative with respect to the root", f_x);
  if (f_x == 0.0) {
    throw_domain_error(function, "partial derivative with respect to the root",
                       f_x, "",
                       " is zero, so the implicit derivative is undefined");
  }

  var result = root;
  reverse_pass_callback([f, result, arena_args_tuple, f_x, msgs]() mutable {
    nested_rev_autodiff nested;
    var residual = math::apply(
        [&](const auto&... arena_args) {
          return eval(f(result.val(), msgs, arena_args...));
        },
        *arena_args_tuple);
    residual.adj() = -result.adj() / f_x;
    grad();
  });

  return result;
}

/**
 * Reverse-mode ITP solver with the default controls.
 */
template <typename F, typename... Args, require_any_st_var<Args...>* = nullptr>
inline var solve_itp(const F& f, double lower, double upper, std::ostream* msgs,
                     const Args&... args) {
  static constexpr const char* function = "solve_itp";
  check_finite(function, "lower bound", lower);
  check_finite(function, "upper bound", upper);
  check_less(function, "lower bound", lower, upper);
  const double width = upper - lower;
  check_finite(function, "upper bound minus lower bound", width);

  return solve_itp_tol(f, lower, upper, 1e-10,
                       internal::itp_default_kappa1(width), 2.0, 1, 200, msgs,
                       args...);
}

}  // namespace math
}  // namespace stan

#endif
