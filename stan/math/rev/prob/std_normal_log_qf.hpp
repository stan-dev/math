#ifndef STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/functor/apply_scalar_ternary.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/prob/std_normal_log_qf.hpp>
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
  return make_callback_var(
      std_normal_log_qf(log_p.val()), [log_p](auto& vi) mutable {
        auto deriv
          = apply_scalar_ternary(
              [](const auto& logp_val, const auto& vi_val, const auto& vi_adj) {
                return vi_adj * exp(logp_val - std_normal_lpdf(vi_val));
              }, log_p.val(), vi.val(), vi.adj());

        as_array_or_scalar(log_p).adj() += as_array_or_scalar(deriv);
      });
}

}  // namespace math
}  // namespace stan
#endif
