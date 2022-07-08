#ifndef STAN_MATH_REV_FUN_OWENS_T_HPP
#define STAN_MATH_REV_FUN_OWENS_T_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam VarVal1 A scalar or Eigen type.
 * @tparam VarVal2 A scalar or Eigen type.
 * @param h var parameter.
 * @param a var parameter.
 * @return The Owen's T function.
 */
template <typename VarVal1, typename VarVal2>
inline auto owens_t(const var_value<VarVal1>& h, const var_value<VarVal2>& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  return make_callback_var(owens_t(h.val(), a.val()), [h_arena, a_arena](auto&& vi) {
    const auto& h_val = as_array_or_scalar(h_arena).val();
    const auto& a_val = as_array_or_scalar(a_arena).val();
    const auto neg_avi_sq_div_2 = -square(h_val) * 0.5;
    const auto one_p_bvi_sq = 1.0 + square(a_val);
    h_arena.adj() += possibly_sum<is_stan_scalar<VarVal1>>(as_array_or_scalar(vi.adj()) * erf(a_val * h_val * INV_SQRT_TWO)
                  * exp(neg_avi_sq_div_2) * INV_SQRT_TWO_PI * -0.5);
    a_arena.adj() += possibly_sum<is_stan_scalar<VarVal2>>(as_array_or_scalar(vi.adj()) * exp(neg_avi_sq_div_2 * one_p_bvi_sq)
                  / (one_p_bvi_sq * TWO_PI));

  });
}

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam VarValue A scalar or Eigen type.
 * @tparam Arith A scalar or Eigen type with an inner arirthmetic scalar value.
 * @param h var parameter.
 * @param a double parameter.
 * @return The Owen's T function.
 */
template <typename VarValue, typename Arith, require_vt_arithmetic<Arith>* = nullptr,
 require_not_std_vector_t<Arith>* = nullptr>
inline auto owens_t(const var_value<VarValue>& h, const Arith& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  using return_type = return_var_matrix_t<decltype(owens_t(h_arena.val(), a_arena)), VarValue, Arith>;
  arena_t<return_type> vi = owens_t(h_arena.val(), a_arena);
  reverse_pass_callback([h_arena, a_arena, vi]() {
    h_arena.adj() += possibly_sum<is_stan_scalar<VarValue>>(as_array_or_scalar(vi.adj()) *
    erf(as_array_or_scalar(a_arena) * as_array_or_scalar(h_arena).val() * INV_SQRT_TWO)
                  * exp(-square(as_array_or_scalar(h_arena).val()) * 0.5) * INV_SQRT_TWO_PI
                  * -0.5);

  });
  return return_type(vi);
}

/**
 * The Owen's T function of h and a.
 *
 * Used to compute the cumulative density function for the skew normal
 * distribution.
 *
 * @tparam VarValue A scalar or Eigen type.
 * @tparam Arith A scalar or Eigen type with an inner arirthmetic scalar value.
 * @param h double parameter.
 * @param a var parameter.
 * @return The Owen's T function.
 */
template <typename Arith, typename VarValue, require_vt_arithmetic<Arith>* = nullptr,
 require_not_std_vector_t<Arith>* = nullptr>
inline auto owens_t(const Arith& h, const var_value<VarValue>& a) {
  auto h_arena = to_arena(h);
  auto a_arena = to_arena(a);
  using return_type = return_var_matrix_t<decltype(owens_t(h_arena, a_arena.val())), VarValue, Arith>;
  arena_t<return_type> vi = owens_t(h_arena, a_arena.val());
  reverse_pass_callback([h_arena, a_arena, vi]() {
    const auto one_p_bvi_sq = eval(1.0 + square(as_array_or_scalar(a_arena.val())));
    a_arena.adj() += possibly_sum<is_stan_scalar<VarValue>>(as_array_or_scalar(vi.adj()) * exp(-0.5 * square(as_array_or_scalar(h_arena)) * one_p_bvi_sq)
                  / (one_p_bvi_sq * TWO_PI));

  });
  return return_type(vi);
}

}  // namespace math
}  // namespace stan
#endif
