#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDZ_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDZ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/functor/apply_scalar_ternary.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to z.
 *
 * @tparam T scalar types of arguments
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @return partial derivative of the incomplete beta with respect to z
 *
 * @pre a > 0
 * @pre b > 0
 * @pre 0 < z <= 1
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr>
inline auto inc_beta_ddz(T1 a, T2 b, T3 z) {
  return exp((b - 1) * log1m(z) + (a - 1) * log(z) + lgamma(a + b) - lgamma(a)
             - lgamma(b));
}

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to z.
 *
 * This overload for double arguments is evaluated by Boost. Where the density
 * is infinite (z = 0 with a < 1, or z = 1 with b < 1) it returns infinity
 * instead of throwing.
 *
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @return partial derivative of the incomplete beta with respect to z
 */
inline double inc_beta_ddz(double a, double b, double z) {
  using boost::math::ibeta_derivative;
  using overflow_ignore
      = boost::math::policies::policy<boost::math::policies::overflow_error<
          boost::math::policies::ignore_error>>;
  return ibeta_derivative(a, b, z, overflow_ignore());
}

/**
 * Enables the vectorized application of the inc_beta_ddz function, when
 *  any arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @tparam T3 type of third input
 * @param a First input
 * @param b Second input
 * @param z Third input
 * @return Partial derivative of the incomplete beta function with respect
 *  to z, applied to the three inputs.
 */
template <typename T1, typename T2, typename T3,
          require_any_container_t<T1, T2, T3>* = nullptr,
          require_all_not_var_matrix_t<T1, T2, T3>* = nullptr>
inline auto inc_beta_ddz(T1&& a, T2&& b, T3&& z) {
  return apply_scalar_ternary(
      [](auto&& d, auto&& e, auto&& f) {
        return inc_beta_ddz(std::forward<decltype(d)>(d),
                            std::forward<decltype(e)>(e),
                            std::forward<decltype(f)>(f));
      },
      std::forward<T1>(a), std::forward<T2>(b), std::forward<T3>(z));
}

}  // namespace math
}  // namespace stan
#endif
