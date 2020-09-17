#ifndef STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the log sum of exponentials.
 *
 * @tparam T_a type of a
 * @tparam T_b type of b
 * @param a first argument
 * @param b first argument
 * @return log of e^a + e^b
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline var log_sum_exp(const T1& a, const T2& b) {
  var res = log_sum_exp(value_of(a), value_of(b));

  reverse_pass_callback([a, b, res]() mutable {
    if (!is_constant<T1>::value)
      forward_as<var>(a).adj()
          += res.adj() * inv_logit(value_of(a) - value_of(b));

    if (!is_constant<T2>::value)
      forward_as<var>(b).adj()
          += res.adj() * inv_logit(value_of(b) - value_of(a));
  });

  return res;
}

/**
 * Returns the log sum of exponentials.
 *
 * @tparam T Type of input vector or matrix.
 * @param x matrix
 */
template <typename T, require_container_st<is_var, T>* = nullptr>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<T>::reduce(x, [](const auto& v) {
    const auto& v_ref = to_ref(v);
    const auto& v_val = value_of(v_ref);

    var res = log_sum_exp(v_val);
    auto arena_v = to_arena(v_ref);

    reverse_pass_callback([arena_v, res]() mutable {
      arena_v.adj()
          += res.adj() * (arena_v.array().val() - res.val()).exp().matrix();
    });

    return res;
  });
}

}  // namespace math
}  // namespace stan
#endif
