#ifndef STAN_MATH_REV_FUNCTOR_ROOT_FINDER_HPP
#define STAN_MATH_REV_FUNCTOR_ROOT_FINDER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_bounded.hpp>
#include <stan/math/prim/functor/root_finder.hpp>
#include <boost/math/tools/roots.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/**
 * var specialization for root solving using Boost's Halley method
 * @tparam FTuple A tuple holding functors whose signatures all match
 *  `(GuessScalar g, Types&&... Args)`.
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
 * @param f_tuple A tuple of functors to calculate the value and any derivates
 * needed.
 * @param guess An initial guess at the root value
 * @param min The minimum possible value for the result, this is used as an
 * initial lower bracket
 * @param max The maximum possible value for the result, this is used as an
 * initial upper bracket
 * @param digits The desired number of binary digits precision
 * @param max_iter An optional maximum number of iterations to perform. On exit,
 * this is updated to the actual number of iterations performed
 * @param args Parameter pack of arguments to pass the the functors in `f_tuple`
 */
template <
    typename FTuple, typename GuessScalar, typename MinScalar,
    typename MaxScalar, typename... Types,
    require_any_st_var<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr,
    require_all_stan_scalar_t<GuessScalar, MinScalar, MaxScalar>* = nullptr>
auto root_finder_tol(FTuple&& f_tuple, const GuessScalar guess,
                     const MinScalar min, const MaxScalar max, const int digits,
                     std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_positive("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  auto arena_args_tuple = make_chainable_ptr(std::make_tuple(eval(args)...));
  auto args_vals_tuple = apply(
      [&](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      *arena_args_tuple);
  // Solve the system
  double theta_dbl = apply(
      [&f_tuple, guess_val = value_of(guess), min_val = value_of(min),
       max_val = value_of(max)](auto&&... vals) {
        return root_finder(f_tuple, guess_val, min_val, max_val, vals...);
      },
      args_vals_tuple);
  double Jf_x;
  double f_x;
  auto&& f = std::get<0>(f_tuple);
  {
    nested_rev_autodiff nested;
    stan::math::var x_var(theta_dbl);
    stan::math::var fx_var = apply(
        [&x_var, &f](auto&&... args) { return f(x_var, std::move(args)...); },
        std::move(args_vals_tuple));
    f_x = fx_var.val();
    fx_var.grad();
    Jf_x = x_var.adj();
  }

  /*
   * Note: Because we put this on the callback stack, if `f` is a lambda
   * its captures must be in Stan's arena memory or trivially destructable.
   */
  return make_callback_var(theta_dbl,
                           [f, arena_args_tuple, Jf_x](auto& ret) mutable {
                             // Eigen::VectorXd eta =
                             // -Jf_x_T_lu_ptr->solve(ret.adj().eval());
                             double eta = -(ret.adj() / Jf_x);
                             // Contract with Jacobian of f with respect to y
                             // using a nested reverse autodiff pass.
                             {
                               nested_rev_autodiff rev;
                               double ret_val = ret.val();
                               auto x_nrad_ = apply(
                                   [&ret_val, &f](const auto&... args) {
                                     return eval(f(ret_val, args...));
                                   },
                                   *arena_args_tuple);
                               x_nrad_.adj() = eta;
                               grad();
                             }
                           });
}

}  // namespace math
}  // namespace stan
#endif
