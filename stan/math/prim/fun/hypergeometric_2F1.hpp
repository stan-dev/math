#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F1_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/asin.hpp>
#include <stan/math/prim/fun/asinh.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <boost/optional.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Calculate the Gauss Hypergeometric (2F1) function for special-case
 * combinations of parameters which can be calculated in closed-form. For
 * more background (and other possible special-cases), see:
 * https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/
 *
 * The return value is wrapped in a boost::optional<> type so that a void
 * return is possible if no special-case rules are applicable
 *
 * @tparam Ta1 Type of scalar first 'a' argument
 * @tparam Ta2 Type of scalar second 'a' argument
 * @tparam Tb Type of scalar 'b' argument
 * @tparam Tz Type of scalar 'z' argument
 * @param[in] a1 First of 'a' arguments to function
 * @param[in] a2 Second of 'a' arguments to function
 * @param[in] b 'b' argument to function
 * @param[in] z Scalar z argument
 * @return Gauss hypergeometric function
 */
template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          typename RtnT = boost::optional<return_type_t<Ta1, Ta1, Tb, Tz>>,
          require_all_arithmetic_t<Ta1, Ta2, Tb, Tz>* = nullptr>
inline RtnT hyper_2F1_special_cases(const Ta1& a1, const Ta2& a2, const Tb& b,
                                    const Tz& z) {
  // https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/01/
  // // NOLINT
  if (z == 0.0) {
    return 1.0;
  }

  // https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/06/01/0001/
  // // NOLINT
  if (a1 == b) {
    return inv(pow(1.0 - z, a2));
  }

  // https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/06/01/0003/
  // // NOLINT
  if (b == (a2 - 1.0)) {
    return (pow((1.0 - z), -a1 - 1.0) * (a2 + z * (a1 - a2 + 1.0) - 1.0))
           / (a2 - 1);
  }

  if (a1 == a2) {
    // https://www.wolframalpha.com/input?i=Hypergeometric2F1%281%2C+1%2C+2%2C+-z%29
    // // NOLINT
    if (a1 == 1.0 && b == 2.0 && z < 0) {
      auto pos_z = abs(z);
      return log1p(pos_z) / pos_z;
    }

    if (a1 == 0.5 && b == 1.5 && z < 1.0) {
      auto sqrt_z = sqrt(abs(z));
      auto numerator
          = (z > 0.0)
                // https://www.wolframalpha.com/input?i=Hypergeometric2F1%281%2F2%2C+1%2F2%2C+3%2F2%2C+z%29
                // // NOLINT
                ? asin(sqrt_z)
                // https://www.wolframalpha.com/input?i=Hypergeometric2F1%281%2F2%2C+1%2F2%2C+3%2F2%2C+-z%29
                // // NOLINT
                : asinh(sqrt_z);
      return numerator / sqrt_z;
    }

    // https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/04/03/
    // // NOLINT
    if (b == (a1 + 1) && z == 0.5) {
      return pow(2, a1 - 1) * a1
             * (digamma((a1 + 1) / 2.0) - digamma(a1 / 2.0));
    }
  }

  if (z == 1.0) {
    // https://www.wolframalpha.com/input?i=Hypergeometric2F1%28a1%2C+a2%2C+a1+%2B+a2+%2B+2%2C+1%29
    // // NOLINT
    if (b == (a1 + a2 + 2)) {
      auto log_2f1 = lgamma(b) - (lgamma(a1 + 2) + lgamma(a2 + 2));
      return exp(log_2f1);
      // https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/02/0001/
      // // NOLINT
    } else if (b > (a1 + a2)) {
      auto log_2f1 = (lgamma(b) + lgamma(b - a1 - a2))
                     - (lgamma(b - a1) + lgamma(b - a2));
      return exp(log_2f1);
    }
  }

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%283%2F2%2C+2%2C+3%2C+-z%29
  // // NOLINT
  if (a1 == 1.5 && a2 == 2.0 && b == 3.0 && z < 0.0) {
    auto abs_z = abs(z);
    auto sqrt_1pz = sqrt(1 + abs_z);
    return -4 * (2 * sqrt_1pz + z - 2) / (sqrt_1pz * square(z));
  }

  return {};
}
}  // namespace internal

/**
 * Returns the Gauss hypergeometric function applied to the
 * input arguments:
 * \f$_2F_1(a_1,a_2;b;z)\f$
 *
 * If the input parameters do not meet convergence criteria, then Euler's
 * transformation is applied to resolve this:
 * https://mathworld.wolfram.com/EulerTransform.html
 *
 * For some special-case combinations of parameters the series is calculated
 * in closed form, see the internal::hyper_2F1_special_cases function for more
 * details.
 *
 * See 'grad_2F1.hpp' for the derivatives wrt each parameter
 *
 * @tparam Ta1 Type of scalar first 'a' argument
 * @tparam Ta2 Type of scalar second 'a' argument
 * @tparam Tb Type of scalar 'b' argument
 * @tparam Tz Type of scalar 'z' argument
 * @param[in] a1 First of 'a' arguments to function
 * @param[in] a2 Second of 'a' arguments to function
 * @param[in] b 'b' argument to function
 * @param[in] z Scalar z argument
 * @return Gauss hypergeometric function
 */
template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          typename ScalarT = return_type_t<Ta1, Ta1, Tb, Tz>,
          typename OptT = boost::optional<ScalarT>,
          require_all_arithmetic_t<Ta1, Ta2, Tb, Tz>* = nullptr>
inline return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1& a1,
                                                          const Ta2& a2,
                                                          const Tb& b,
                                                          const Tz& z) {
  check_finite("hypergeometric_2F1", "a1", a1);
  check_finite("hypergeometric_2F1", "a2", a2);
  check_finite("hypergeometric_2F1", "b", b);
  check_finite("hypergeometric_2F1", "z", z);

  check_not_nan("hypergeometric_2F1", "a1", a1);
  check_not_nan("hypergeometric_2F1", "a2", a2);
  check_not_nan("hypergeometric_2F1", "b", b);
  check_not_nan("hypergeometric_2F1", "z", z);

  // Check whether value can be calculated by any special-case rules
  // before estimating infinite sum
  OptT special_case_a1a2 = internal::hyper_2F1_special_cases(a1, a2, b, z);
  if (special_case_a1a2.is_initialized()) {
    return special_case_a1a2.get();
  }

  // Check whether any special case rules apply with 'a' arguments reversed
  // as 2F1(a1, a2, b, z) = 2F1(a2, a1, b, z)
  OptT special_case_a2a1 = internal::hyper_2F1_special_cases(a2, a1, b, z);
  if (special_case_a2a1.is_initialized()) {
    return special_case_a2a1.get();
  }

  Eigen::Matrix<double, 2, 1> a_args(2);
  Eigen::Matrix<double, 1, 1> b_args(1);

  try {
    check_2F1_converges("hypergeometric_2F1", a1, a2, b, z);

    a_args << a1, a2;
    b_args << b;
    return hypergeometric_pFq(a_args, b_args, z);
  } catch (const std::exception& e) {
    // Apply Euler's hypergeometric transformation if function
    // will not converge with current arguments
    ScalarT a1_t = b - a1;
    ScalarT a2_t = a2;
    ScalarT b_t = b;
    ScalarT z_t = z / (z - 1);

    check_2F1_converges("hypergeometric_2F1", a1_t, a2_t, b_t, z_t);

    a_args << a1_t, a2_t;
    b_args << b_t;
    return hypergeometric_pFq(a_args, b_args, z_t) / pow(1 - z, a2);
  }
}
}  // namespace math
}  // namespace stan
#endif
