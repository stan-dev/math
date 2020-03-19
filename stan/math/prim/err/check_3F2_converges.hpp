#ifndef STAN_MATH_PRIM_ERR_CHECK_3F2_CONVERGES_HPP
#define STAN_MATH_PRIM_ERR_CHECK_3F2_CONVERGES_HPP

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
 * Check if the hypergeometric function (3F2) called with
 * supplied arguments will converge, assuming arguments are
 * finite values.
 * @tparam T_a1 Type of a1
 * @tparam T_a2 Type of a2
 * @tparam T_a3 Type of a3
 * @tparam T_b1 Type of b1
 * @tparam T_b2 Type of b2
 * @tparam T_z Type of z
 * @param function Name of function ultimately relying on 3F2 (for error
 &   messages)
 * @param a1 Variable to check
 * @param a2 Variable to check
 * @param a3 Variable to check
 * @param b1 Variable to check
 * @param b2 Variable to check
 * @param z Variable to check
 * @throw <code>domain_error</code> if 3F2(a1, a2, a3, b1, b2, z)
 *   does not meet convergence conditions, or if any coefficient is NaN.
 */
template <typename T_a1, typename T_a2, typename T_a3, typename T_b1,
          typename T_b2, typename T_z>
inline void check_3F2_converges(const char* function, const T_a1& a1,
                                const T_a2& a2, const T_a3& a3, const T_b1& b1,
                                const T_b2& b2, const T_z& z) {
  using std::fabs;
  using std::floor;

  double a1d = value_of_rec(a1);
  double a2d = value_of_rec(a2);
  double a3d = value_of_rec(a3);
  double b1d = value_of_rec(b1);
  double b2d = value_of_rec(b2);
  double zd = value_of_rec(z);

  check_not_nan("check_3F2_converges", "a1", a1d);
  check_not_nan("check_3F2_converges", "a2", a2d);
  check_not_nan("check_3F2_converges", "a3", a3d);
  check_not_nan("check_3F2_converges", "b1", b1d);
  check_not_nan("check_3F2_converges", "b2", b2d);
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
  double fabs_a2d = fabs(a3d);
  if (is_nonpositive_integer(a3d) && fabs_a3d >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs_a3d);
  }

  bool is_undefined
      = (is_nonpositive_integer(b1d) && fabs(b1d) <= num_terms)
        || (is_nonpositive_integer(b2d) && fabs(b2d) <= num_terms);

  if (is_polynomial && !is_undefined) {
    return;
  }
  double fabs_zd = fabs(zd);
  if (fabs_zd < 1.0 && !is_undefined) {
    return;
  }
  if (fabs_zd == 1.0 && !is_undefined && b1d + b2d > a1d + a2d + a3d) {
    return;
  }

  std::stringstream msg;
  msg << "called from function '" << function << "', "
      << "hypergeometric function 3F2 does not meet convergence "
      << "conditions with given arguments. "
      << "a1: " << a1d << ", a2: " << a2d << ", a3: " << a3d << ", b1: " << b1d
      << ", b2: " << b2d << ", z: " << zd;
  throw std::domain_error(msg.str());
}

}  // namespace math
}  // namespace stan
#endif
