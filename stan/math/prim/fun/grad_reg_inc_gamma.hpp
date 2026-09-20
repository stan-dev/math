#ifndef STAN_MATH_PRIM_FUN_GRAD_REG_INC_GAMMA_HPP
#define STAN_MATH_PRIM_FUN_GRAD_REG_INC_GAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <unsupported/Eigen/SpecialFunctions>
#include <limits>

namespace stan {
namespace math {

/**
 * Gradient of the regularized incomplete gamma function igamma(a, z) with
 * respect to the shape parameter, d/da Q(a, z), where Q is the upper
 * regularized incomplete gamma function.
 *
 * Delegates to `Eigen::numext::igamma_der_a`, which returns d/da of the
 * lower P(a, z); the sign is flipped here. For autodiff scalar types the
 * Cephes helpers that Eigen needs are specialized in
 * `stan/math/fwd/fun/Eigen_SpecialFunctions.hpp` and
 * `stan/math/rev/fun/Eigen_SpecialFunctions.hpp`.
 *
 * `g`, `dig`, `precision` and `max_steps` are accepted for signature
 * compatibility and are not used.
 *
 * @tparam T1 type of the shape parameter
 * @tparam T2 type of the location parameter
 * @param a shape parameter, a > 0
 * @param z location z >= 0
 * @param g unused
 * @param dig unused
 * @param precision unused
 * @param max_steps unused
 * @return d/da of the upper regularized incomplete gamma function
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> grad_reg_inc_gamma(T1 a, T2 z, T1 g, T1 dig,
                                                double precision = 1e-6,
                                                int max_steps = 1e5) {
  using TP = return_type_t<T1, T2>;

  if (is_any_nan(a, z, g, dig)) {
    return std::numeric_limits<TP>::quiet_NaN();
  }

  return -Eigen::numext::igamma_der_a(TP(a), TP(z));
}

}  // namespace math
}  // namespace stan
#endif
