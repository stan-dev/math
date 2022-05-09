#ifndef STAN_MATH_REV_FUN_INV_PHI_LOG_HPP
#define STAN_MATH_REV_FUN_INV_PHI_LOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv_Phi_log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse of unit normal cumulative density function.
 *
 * The derivative is the reciprocal of unit normal density function,
 *
 * @param log_p log probability
 * @return the unit normal inverse cdf evaluated at log_p
 */
template <typename T, require_all_stan_scalar_t<T>* = nullptr,
          require_any_var_t<T>* = nullptr>
inline var inv_Phi_log(const T& log_p) {
  double log_p_val = value_of(log_p);
  double z = inv_Phi_log(log_p_val);

  return make_callback_var(z, [log_p, log_p_val](auto& vi) mutable {
    log_p.adj() += log_p_val * vi.adj() * SQRT_TWO_PI
                       / std::exp(-0.5 * vi.val() * vi.val())
                   + vi.val() / exp(log_p_val);
  });
}

/**
 * Return the elementwise inverse of unit normal cumulative density function.
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param p Probability vector
 * @return Elementwise unit normal inverse cdf
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto inv_Phi_log(const T& log_p) {
  auto log_p_val = value_of(log_p);
  auto z = inv_Phi_log(log_p_val);
  return make_callback_var(
      inv_Phi_log(log_p.val()), [log_p, log_p_val](auto& vi) mutable {
        as_array_or_scalar(log_p.adj())
            += log_p_val.array() * vi.val().array().exp() * vi.adj().array()
                   * SQRT_TWO_PI / (-0.5 * vi.val().array().square()).exp()
               + vi.val().array() / log_p_val.array().exp();
      });
}

}  // namespace math
}  // namespace stan
#endif
