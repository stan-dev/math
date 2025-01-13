#ifndef STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/prob/std_normal_log_qf.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <cmath>

namespace stan {
namespace math {
/**
 * Return the elementwise inverse of unit normal cumulative density function.
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param log_p log probability vector
 * @return Elementwise unit normal inverse cdf
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto std_normal_log_qf(const var_value<T>& log_p) {
  const auto& arena_rtn = to_arena(std_normal_log_qf(log_p.val()));
  return make_callback_var(arena_rtn, [log_p, arena_rtn](auto& vi) mutable {
    auto deriv = apply_scalar_binary(
        log_p.val(), arena_rtn, [](const auto& logp_val, const auto& rtn_val) {
          return exp(logp_val - std_normal_lpdf(rtn_val));
        });
    log_p.adj() += elt_multiply(vi.adj(), deriv);
  });
}

}  // namespace math
}  // namespace stan
#endif
