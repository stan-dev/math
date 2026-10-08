#ifndef STAN_MATH_PRIM_FUN_INC_BETA_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/apply_scalar_ternary.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/beta.hpp>

namespace stan {
namespace math {

/**
 * The normalized incomplete beta function of a, b, with outcome x.
 *
 * Used to compute the cumulative density function for the beta
 * distribution.
 *
 * @param a Shape parameter a >= 0; a and b can't both be 0
 * @param b Shape parameter b >= 0
 * @param x Random variate. 0 <= x <= 1
 * @throws if constraints are violated or if any argument is NaN
 * @return The normalized incomplete beta function.
 */
inline double inc_beta(double a, double b, double x) {
  check_not_nan("inc_beta", "a", a);
  check_not_nan("inc_beta", "b", b);
  check_not_nan("inc_beta", "x", x);
  return boost::math::ibeta(a, b, x, boost_policy_t<>());
}

namespace internal {

/**
 * Return the complement of the regularized incomplete beta function,
 * 1 - I_x(a, b) = I_{1-x}(b, a), from both x and 1 - x.
 *
 * Boost forms 1 - x from its argument, so an argument close to 1 loses
 * relative precision in 1 - x. This overload calls Boost with whichever
 * of x and 1 - x is at most 1/2: ibetac(a, b, x) or ibeta(b, a, 1 - x).
 * Both keep their relative precision when the complement is below eps.
 *
 * Below the mean a / (a + b), while I_x(a, b) is at most 1/2, the
 * complement is formed as 1 - I_x(a, b), which then has full precision.
 * This also avoids Boost's arcsine case a = b = 1/2, where ibetac
 * evaluates asin(sqrt(1 - x)) and loses digits for small x.
 *
 * @param a first shape, a > 0
 * @param b second shape, b > 0
 * @param x argument, 0 <= x <= 1
 * @param one_m_x 1 - x, computed by the caller without cancellation
 * @return 1 - I_x(a, b)
 */
inline double inc_beta_complement(double a, double b, double x,
                                  double one_m_x) {
  check_not_nan("inc_beta_complement", "a", a);
  check_not_nan("inc_beta_complement", "b", b);
  check_not_nan("inc_beta_complement", "x", x);
  if (x > 0.5) {
    return boost::math::ibeta(b, a, one_m_x, boost_policy_t<>());
  }
  if (x < a / (a + b)) {
    const double inc = boost::math::ibeta(a, b, x, boost_policy_t<>());
    if (inc <= 0.5) {
      return 1.0 - inc;
    }
  }
  return boost::math::ibetac(a, b, x, boost_policy_t<>());
}

/**
 * Return the complement of the regularized incomplete beta function for
 * autodiff arguments, by the symmetry relation through inc_beta.
 *
 * @tparam T autodiff type
 * @param a first shape, a > 0
 * @param b second shape, b > 0
 * @param x argument, 0 <= x <= 1; not used
 * @param one_m_x 1 - x
 * @return 1 - I_x(a, b)
 */
template <typename T>
inline T inc_beta_complement(const T& a, const T& b, const T& x,
                             const T& one_m_x) {
  return inc_beta(b, a, one_m_x);
}

}  // namespace internal

/**
 * Enables the vectorized application of the inc_beta function, when
 *  any arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @tparam T3 type of third input
 * @param a First input
 * @param b Second input
 * @param c Third input
 * @return Incomplete Beta function applied to the three inputs.
 */
template <typename T1, typename T2, typename T3,
          require_any_container_t<T1, T2, T3>* = nullptr,
          require_all_not_var_matrix_t<T1, T2, T3>* = nullptr>
inline auto inc_beta(T1&& a, T2&& b, T3&& c) {
  return apply_scalar_ternary(
      [](auto&& d, auto&& e, auto&& f) {
        return inc_beta(std::forward<decltype(d)>(d),
                        std::forward<decltype(e)>(e),
                        std::forward<decltype(f)>(f));
      },
      std::forward<T1>(a), std::forward<T2>(b), std::forward<T3>(c));
}

}  // namespace math
}  // namespace stan
#endif
