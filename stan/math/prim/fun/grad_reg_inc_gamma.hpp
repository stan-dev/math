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
 * The work is done by `Eigen::numext::igamma_der_a`, which differentiates
 * the Cephes power series and the Cephes continued fraction term by term
 * and stops at machine epsilon. It returns d/da of the lower P(a, z), so
 * the sign is flipped here.
 *
 * Eigen's implementation is generic in the scalar type. Stan supplies the
 * three Cephes helpers it needs for autodiff scalars in
 * `stan/math/fwd/fun/Eigen_SpecialFunctions.hpp` and
 * `stan/math/rev/core/Eigen_SpecialFunctions.hpp`, so every autodiff order
 * uses this same algorithm.
 *
 * `g`, `dig`, `precision` and `max_steps` are accepted and ignored. They
 * belonged to the previous hand-written series. The signature is unchanged
 * so that the 27 call sites do not change; removing the arguments is
 * proposed separately.
 *
 * @tparam T1 type of the shape parameter
 * @tparam T2 type of the location parameter
 * @param a shape parameter, a > 0
 * @param z location z >= 0
 * @param g ignored; previously stan::math::tgamma(a)
 * @param dig ignored; previously stan::math::digamma(a)
 * @param precision ignored; previously the series tolerance
 * @param max_steps ignored; previously the series iteration limit
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
