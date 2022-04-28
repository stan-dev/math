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

/*
// TODO: Can just have one signature with another function for constructing
// tuple
template <
  typename F, typename FDiv, typename FDivDiv, typename GuessScalar,
  typename MinScalar, typename MaxScalar, typename... Types,
  require_any_st_var<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr>
auto root_finder_tol(F&& f, FDiv&& f_div, FDivDiv&& f_div_div,
                   GuessScalar guess, MinScalar min, MaxScalar max,
                   int digits, std::uintmax_t& max_iter, Types&&... args) {
check_bounded("root_finder", "initial guess", guess, min, max);
return_type_t<GuessScalar> ret = 0;
auto f_plus_div = [&f, &f_div, &f_div_div, &args...](auto&& g) {
  return std::make_tuple(f(g, args...), f_div(g, args...),
                         f_div_div(g, args...));
};
try {
  ret = boost::math::tools::halley_iterate(f_plus_div, guess, min, max,
                                           digits, max_iter);
} catch (const std::exception& e) {
  throw e;
}
return ret;
}
*/

template <
    typename F, typename FDiv, typename FDivDiv, typename GuessScalar,
    typename MinScalar, typename MaxScalar, typename... Types,
    require_any_st_var<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr>
auto root_finder_tol(F&& f, FDiv&& f_div, FDivDiv&& f_div_div,
                     GuessScalar guess, MinScalar min, MaxScalar max,
                     int digits, std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_nonnegative("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  auto arena_args_tuple = make_chainable_ptr(std::make_tuple(eval(args)...));
  auto args_vals_tuple = apply(
      [&](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      *arena_args_tuple);
  // Solve the system
  double theta_dbl = apply(
      [&](auto&&... vals) {
        return root_finder(f, f_div, f_div_div, value_of(guess), value_of(min),
                           value_of(max), vals...);
      },
      args_vals_tuple);
  double Jf_x;
  double f_x;
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
  stan::math::var ret = theta_dbl;
  reverse_pass_callback([f, ret, arena_args_tuple, Jf_x]() mutable {
    // Eigen::VectorXd eta = -Jf_x_T_lu_ptr->solve(ret.adj().eval());
    double eta = -(ret.adj() / Jf_x);
    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
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

  return ret;
}

}  // namespace math
}  // namespace stan
#endif
