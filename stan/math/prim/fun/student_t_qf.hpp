#ifndef STAN_MATH_PRIM_FUN_STUDENT_T_QF_HPP
#define STAN_MATH_PRIM_FUN_STUDENT_T_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/detail/t_distribution_inv.hpp>
#include <boost/math/distributions/students_t.hpp>

namespace stan {
namespace math {

/**
 * The quantile function of the Student's T density function for the
 * given degrees of freedom and probability value.
 *
 * @param p random variate. 0 <= p <= 1
 * @param df degrees of freedom parameter df > 0
 * @throws if constraints are violated or if any argument is NaN
 * @return real value of the inverse cdf for the Student's T distribution.
 */
inline double student_t_qf(double p, double df) {
  check_not_nan("student_t_qf", "df", df);
  check_not_nan("student_t_qf", "p", p);
  check_positive("student_t_qf", "df", df);
  check_bounded("student_t_qf", "p", p, 0, 1);
  if (p == 0) {
      return stan::math::NEGATIVE_INFTY;
  }
  if (p == 1) {
      return stan::math::INFTY;
  }
  if (p == 0.5) {
      return 0;
  }
  boost::math::students_t dist(df);
  return boost::math::quantile(dist, p);
}

/**
 * Enables the vectorized version of student_t_qf() when the 
 * first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Student's T quantile function applied to the two inputs
 */
template <typename T1, typename T2,
    require_any_container_t<T1, T2>* = nullptr>
inline auto student_t_qf(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return student_t_qf(c, d);
  }); 
  }

}  // namespace math
}  // namespace stan
#endif
