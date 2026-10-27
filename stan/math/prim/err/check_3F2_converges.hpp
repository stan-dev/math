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

  check_not_nan("check_3F2_converges", "a1", a1);
  check_not_nan("check_3F2_converges", "a2", a2);
  check_not_nan("check_3F2_converges", "a3", a3);
  check_not_nan("check_3F2_converges", "b1", b1);
  check_not_nan("check_3F2_converges", "b2", b2);
  check_not_nan("check_3F2_converges", "z", z);

  // A numerator parameter a = -m (m a nonnegative integer) makes the terms
  // zero from term m + 1 on, so the series is a polynomial that ends at the
  // smallest such m.
  bool is_polynomial = false;
  double num_terms = 0;
  auto add_end = [&](const auto& a) {
    if (is_nonpositive_integer(a)) {
      const double m = fabs(value_of_rec(a));
      num_terms = is_polynomial ? std::fmin(num_terms, m) : m;
      is_polynomial = true;
    }
  };
  add_end(a1);
  add_end(a2);
  add_end(a3);

  // A denominator parameter b = -m (m a nonnegative integer) is a pole from
  // term m + 1 on, because (b)_k is zero for k > m. A polynomial ends at
  // term num_terms, so there b is a pole only if m < num_terms; in an
  // infinite series b is always a pole.
  auto is_pole = [&](const auto& b) {
    return is_nonpositive_integer(b)
           && (!is_polynomial || fabs(value_of_rec(b)) < num_terms);
  };
  bool is_undefined = is_pole(b1) || is_pole(b2);

  if (is_polynomial && !is_undefined) {
    return;
  }
  if (fabs(z) < 1.0 && !is_undefined) {
    return;
  }
  if (fabs(z) == 1.0 && !is_undefined && b1 + b2 > a1 + a2 + a3) {
    return;
  }

  std::stringstream msg;
  msg << "called from function '" << function << "', "
      << "hypergeometric function 3F2 does not meet convergence "
      << "conditions with given arguments. "
      << "a1: " << a1 << ", a2: " << a2 << ", a3: " << a3 << ", b1: " << b1
      << ", b2: " << b2 << ", z: " << z;
  throw std::domain_error(msg.str());
}

}  // namespace math
}  // namespace stan
#endif
