#ifndef STAN_MATH_PRIM_FUN_LOWER_REG_INC_GAMMA_HPP
#define STAN_MATH_PRIM_FUN_LOWER_REG_INC_GAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <unsupported/Eigen/SpecialFunctions>
#include <limits>

namespace stan {
namespace math {

/**
 * Computes the gradient of the lower regularized incomplete gamma
 * function, d/da P(a, z).
 *
 * Delegates to `Eigen::numext::igamma_der_a`, which returns exactly this
 * quantity. For autodiff scalar types the Cephes helpers that Eigen needs
 * are specialized in `stan/math/fwd/fun/Eigen_SpecialFunctions.hpp` and
 * `stan/math/rev/fun/Eigen_SpecialFunctions.hpp`.
 *
 * `precision` and `max_steps` are accepted for signature compatibility and
 * are not used.
 *
 * Infinite z throws, to match the behaviour of, for example,
 * boost::math::gamma_p.
 *
 * @tparam T1 type of a
 * @tparam T2 type of z
 * @param[in] a shared with complete Gamma, a > 0
 * @param[in] z value to integrate up to, z >= 0
 * @param[in] precision unused
 * @param[in] max_steps unused
 * @return d/da of the lower regularized incomplete gamma function
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> grad_reg_lower_inc_gamma(const T1& a, const T2& z,
                                                      double precision = 1e-10,
                                                      int max_steps = 1e5) {
  using TP = return_type_t<T1, T2>;

  if (is_any_nan(a, z)) {
    return std::numeric_limits<TP>::quiet_NaN();
  }

  check_positive_finite("grad_reg_lower_inc_gamma", "a", a);

  if (z == 0.0) {
    return 0.0;
  }
  check_positive_finite("grad_reg_lower_inc_gamma", "z", z);

  return Eigen::numext::igamma_der_a(TP(a), TP(z));
}

}  // namespace math
}  // namespace stan

#endif
