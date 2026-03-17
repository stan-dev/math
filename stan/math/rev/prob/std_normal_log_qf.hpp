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
template <typename T, require_var_t<T>* = nullptr>
inline auto std_normal_log_qf(T&& log_p) {
  auto arena_rtn = to_arena(std_normal_log_qf(log_p.val()));
  return make_callback_var(arena_rtn, [log_p, arena_rtn](auto& vi) mutable {
    if constexpr (is_eigen<decltype(arena_rtn)>::value) {
      auto deriv = exp(log_p.val() - arena_rtn.unaryExpr([](auto x) {
        return std_normal_lpdf(x);
      }));
      log_p.adj() += elt_multiply(vi.adj(), deriv);
    } else {
      auto deriv = exp(log_p.val() - std_normal_lpdf(arena_rtn));
      log_p.adj() += vi.adj() * deriv;
    }
  });
}

}  // namespace math
}  // namespace stan
#endif
