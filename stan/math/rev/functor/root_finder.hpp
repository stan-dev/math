#ifndef STAN_MATH_REV_FUNCTOR_ROOT_FINDER_HPP
#define STAN_MATH_REV_FUNCTOR_ROOT_FINDER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_bounded.hpp>
#include <stan/math/prim/functor/root_finder.hpp>
#include <stan/math/rev/core/nested_rev_autodiff.hpp>
#include <stan/math/rev/meta.hpp>
#include <boost/math/tools/roots.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/**
 * var specialization for root solving using Boost's Halley method
 * @tparam FRootFunc A struct or class with a static function called `run`.
 *  The structs `run` function must have a boolean template parameter that
 *  when `true` returns a tuple containing the function result and the
 * derivatives needed for the root finder. When the boolean template parameter
 * is `false` the function should return a single value containing the function
 * result.
 * @tparam SolverFun One of the three struct types used to call the root solver.
 *  (`NewtonRootSolver`, `HalleyRootSolver`, `SchroderRootSolver`).
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
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
    typename FRootFunc, typename SolverFun, typename GuessScalar,
    typename MinScalar, typename MaxScalar, typename... Types,
    require_any_st_var<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr,
    require_all_stan_scalar_t<GuessScalar, MinScalar, MaxScalar>* = nullptr,
    require_any_not_stan_scalar_t<Types...>* = nullptr>
inline auto root_finder_tol(const GuessScalar guess, const MinScalar min,
                            const MaxScalar max, const int digits,
                            std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_positive("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  auto arena_args_tuple
      = make_chainable_ptr(std::make_tuple(eval(std::forward<Types>(args))...));
  auto args_vals_tuple = apply(
      [&](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      *arena_args_tuple);
  // Solve the system
  double theta_dbl = apply(
      [&max_iter, digits, guess_val = value_of(guess), min_val = value_of(min),
       max_val = value_of(max)](auto&&... vals) {
        return root_finder_tol<FRootFunc, SolverFun>(
            guess_val, min_val, max_val, digits > 20 ? digits : 21, max_iter,
            vals...);
      },
      args_vals_tuple);
  double Jf_x;
  {
    nested_rev_autodiff nested;
    stan::math::var x_var(theta_dbl);
    stan::math::var fx_var = apply(
        [&x_var](auto&&... args) {
          return std::decay_t<FRootFunc>::template run<false>(
              x_var, std::move(args)...);
        },
        std::move(args_vals_tuple));
    fx_var.grad();
    Jf_x = x_var.adj();
  }

  /*
   * Note: Because we put this on the callback stack, if `f` is a lambda
   * its captures must be in Stan's arena memory or trivially destructable.
   */
  return make_callback_var(
      theta_dbl, [arena_args_tuple, Jf_x](auto& ret) mutable {
        {
          nested_rev_autodiff rev;
          double eta = -(ret.adj() / Jf_x);
          double ret_val = ret.val();
          auto x_nrad_ = apply(
              [&ret_val](const auto&... args) {
                auto f = internal::make_root_func<false, FRootFunc>();
                return eval(std::decay_t<FRootFunc>::template run<false>(
                    ret_val, args...));
              },
              *arena_args_tuple);
          x_nrad_.adj() = eta;
          grad();
        }
      });
}

/**
 * var specialization for root solving using Boost's Halley method
 * @tparam FRootFunc A struct or class with a static function called `run`.
 *  The structs `run` function must have a boolean template parameter that
 *  when `true` returns a tuple containing the function result and the
 * derivatives needed for the root finder. When the boolean template parameter
 * is `false` the function should return a single value containing the function
 * result.
 * @tparam SolverFun One of the three struct types used to call the root solver.
 *  (`NewtonRootSolver`, `HalleyRootSolver`, `SchroderRootSolver`).
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
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
    typename FRootFunc, typename SolverFun, typename GuessScalar,
    typename MinScalar, typename MaxScalar, typename... Types,
    require_any_st_var<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr,
    require_all_stan_scalar_t<GuessScalar, MinScalar, MaxScalar>* = nullptr,
    require_all_stan_scalar_t<Types...>* = nullptr>
inline auto root_finder_tol(const GuessScalar guess, const MinScalar min,
                            const MaxScalar max, const int digits,
                            std::uintmax_t& max_iter, Types... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_positive("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  // Solve the system
  double theta_dbl = root_finder_tol<FRootFunc, SolverFun>(
      value_of(guess), value_of(min), value_of(max), digits > 20 ? digits : 21,
      max_iter, value_of(args)...);
  double Jf_x;
  {
    nested_rev_autodiff nested;
    stan::math::var x_var(theta_dbl);
    stan::math::var fx_var = std::decay_t<FRootFunc>::template run<false>(
        x_var, value_of(args)...);
    fx_var.grad();
    Jf_x = x_var.adj();
  }

  /*
   * Note: Because we put this on the callback stack, if `f` is a lambda
   * its captures must be in Stan's arena memory or trivially destructable.
   */
  return make_callback_var(theta_dbl, [Jf_x, args...](auto& ret) mutable {
    {
      nested_rev_autodiff rev;
      double eta = -(ret.adj() / Jf_x);
      double ret_val = ret.val();
      auto f = internal::make_root_func<false, FRootFunc>();
      auto x_nrad
          = std::decay_t<FRootFunc>::template run<false>(ret_val, args...);
      x_nrad.adj() = eta;
      grad();
    }
  });
}

}  // namespace math
}  // namespace stan
#endif
