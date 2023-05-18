#ifndef STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_REV_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/sign.hpp>
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
        auto vi_array = as_array_or_scalar(vi.val());
        auto vi_sign = sign(as_array_or_scalar(vi.adj()));

        const auto& deriv = as_array_or_scalar(log_p).val()
                            + log(as_array_or_scalar(vi.adj()) * vi_sign)
                            - NEG_LOG_SQRT_TWO_PI + 0.5 * square(vi_array);
        as_array_or_scalar(log_p).adj() += vi_sign * exp(deriv);
      });
}

}  // namespace math
}  // namespace stan
#endif
