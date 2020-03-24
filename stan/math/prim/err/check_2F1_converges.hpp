#ifndef STAN_MATH_PRIM_ERR_CHECK_2F1_CONVERGES_HPP
#define STAN_MATH_PRIM_ERR_CHECK_2F1_CONVERGES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <stan/math/prim/fun/is_nonpositive_integer.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Check if the hypergeometric function (2F1) called with
 * supplied arguments will converge, assuming arguments are
 * finite values.
 * @tparam T_a1 Type of a1
 * @tparam T_a2 Type of a2
 * @tparam T_b1 Type of b1
 * @tparam T_z Type of z
 * @param function Name of function ultimately relying on 2F1 (for error
 *   messages)
 * @param a1 Variable to check
 * @param a2 Variable to check
 * @param b1 Variable to check
 * @param z Variable to check
 * @throw <code>domain_error</code> if 2F1(a1, a2, b1, z)
 *   does not meet convergence conditions, or if any coefficient is NaN.
 */
template <typename T_a1, typename T_a2, typename T_b1, typename T_z>
inline void check_2F1_converges(const char* function, const T_a1& a1,
                                const T_a2& a2, const T_b1& b1, const T_z& z) {
  using std::fabs;
  using std::floor;

  double a1d = value_of_rec(a1);
  double a2d = value_of_rec(a2);
  double b1d = value_of_rec(b1);
  double zd = value_of_rec(z);

  check_not_nan("check_3F2_converges", "a1", a1d);
  check_not_nan("check_3F2_converges", "a2", a2d);
  check_not_nan("check_3F2_converges", "b1", b1d);
  check_not_nan("check_3F2_converges", "z", zd);

  int num_terms = 0;
  bool is_polynomial = false;

  double fabs_a1d = fabs(a1d);
  if (is_nonpositive_integer(a1d) && fabs_a1d >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs_a1d);
  }
  double fabs_a2d = fabs(a2d);
  if (is_nonpositive_integer(a2d) && fabs_a2d >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs_a2d);
  }

  bool is_undefined = is_nonpositive_integer(b1d) && fabs(b1d) <= num_terms;

  double fabs_zd = fabs(zd);
  if (!is_undefined
      && (is_polynomial || fabs_zd < 1 || (fabs_zd == 1 && b1d > a1d + a2d))) {
    return;
  }

  std::stringstream msg;
  msg << "called from function '" << function << "', "
      << "hypergeometric function 2F1 does not meet convergence "
      << "conditions with given arguments. "
      << "a1: " << a1d << ", a2: " << a2d << ", "
      << "b1: " << b1d << ", z: " << zd;
  throw std::domain_error(msg.str());
}

}  // namespace math
}  // namespace stan
#endif
