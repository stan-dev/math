#ifndef STAN_MATH_REV_FUN_OWENS_T_HPP
#define STAN_MATH_REV_FUN_OWENS_T_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/erf.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam Var1 A scalar or Eigen type whose `scalar_type` is an var.
 * @tparam Var2 A scalar or Eigen type whose `scalar_type` is an var.
 * @param h var parameter.
 * @param a var parameter.
 * @return The Owen's T function.
 */
template <typename Var1, typename Var2,
          require_all_st_var<Var1, Var2>* = nullptr,
          require_all_not_std_vector_t<Var1, Var2>* = nullptr>
inline auto owens_t(const Var1& h, const Var2& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  using return_type
      = return_var_matrix_t<decltype(owens_t(h_arena.val(), a_arena.val())),
                            Var1, Var2>;
  arena_t<return_type> ret = owens_t(h_arena.val(), a_arena.val());
  reverse_pass_callback([h_arena, a_arena, ret]() mutable {
    const auto& h_val = as_value_array_or_scalar(h_arena);
    const auto& a_val = as_value_array_or_scalar(a_arena);
    const auto neg_h_sq_div_2 = stan::math::eval(-square(h_val) * 0.5);
    const auto one_p_a_sq = stan::math::eval(1.0 + square(a_val));
    as_array_or_scalar(h_arena).adj() += possibly_sum<is_stan_scalar<Var1>>(
        as_array_or_scalar(ret.adj()) * erf(a_val * h_val * INV_SQRT_TWO)
        * exp(neg_h_sq_div_2) * INV_SQRT_TWO_PI * -0.5);
    as_array_or_scalar(a_arena).adj() += possibly_sum<is_stan_scalar<Var2>>(
        as_array_or_scalar(ret.adj()) * exp(neg_h_sq_div_2 * one_p_a_sq)
        / (one_p_a_sq * TWO_PI));
  });
  return return_type(ret);
}

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam Var A scalar or Eigen type whose `scalar_type` is an var.
 * @tparam Arith A scalar or Eigen type with an inner arirthmetic scalar value.
 * @param h var parameter.
 * @param a double parameter.
 * @return The Owen's T function.
 */
template <typename Var, typename Arith, require_st_arithmetic<Arith>* = nullptr,
          require_all_not_std_vector_t<Var, Arith>* = nullptr,
          require_st_var<Var>* = nullptr>
inline auto owens_t(const Var& h, const Arith& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  using return_type
      = return_var_matrix_t<decltype(owens_t(h_arena.val(), a_arena)), Var,
                            Arith>;
  arena_t<return_type> ret = owens_t(h_arena.val(), a_arena);
  reverse_pass_callback([h_arena, a_arena, ret]() mutable {
    const auto& h_val = as_value_array_or_scalar(h_arena);
    as_array_or_scalar(h_arena).adj() += possibly_sum<is_stan_scalar<Var>>(
        as_array_or_scalar(ret.adj())
        * erf(as_array_or_scalar(a_arena) * h_val * INV_SQRT_TWO)
        * exp(-square(h_val) * 0.5) * INV_SQRT_TWO_PI * -0.5);
  });
  return return_type(ret);
}

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam Var A scalar or Eigen type whose `scalar_type` is an var.
 * @tparam Arith A scalar or Eigen type with an inner arithmetic scalar value.
 * @param h double parameter.
 * @param a var parameter.
 * @return The Owen's T function.
 */
template <typename Arith, typename Var, require_st_arithmetic<Arith>* = nullptr,
          require_all_not_std_vector_t<Var, Arith>* = nullptr,
          require_st_var<Var>* = nullptr>
inline auto owens_t(const Arith& h, const Var& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  using return_type
      = return_var_matrix_t<decltype(owens_t(h_arena, a_arena.val())), Var,
                            Arith>;
  arena_t<return_type> ret = owens_t(h_arena, a_arena.val());
  reverse_pass_callback([h_arena, a_arena, ret]() mutable {
    const auto one_p_a_sq
        = eval(1.0 + square(as_value_array_or_scalar(a_arena)));
    as_array_or_scalar(a_arena).adj() += possibly_sum<is_stan_scalar<Var>>(
        as_array_or_scalar(ret.adj())
        * exp(-0.5 * square(as_array_or_scalar(h_arena)) * one_p_a_sq)
        / (one_p_a_sq * TWO_PI));
  });
  return return_type(ret);
}

}  // namespace math
}  // namespace stan
#endif
