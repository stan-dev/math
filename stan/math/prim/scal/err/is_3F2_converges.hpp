#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_3F2_CONVERGES_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_3F2_CONVERGES_HPP

#include <stan/math/prim/scal/err/is_not_nan.hpp>
#include <stan/math/prim/scal/fun/is_nonpositive_integer.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Check if the hypergeometric function (3F2) called with supplied arguments
 * will converge, assuming arguments are finite values.
 * @tparam T_a1 Type of `a1`
 * @tparam T_a2 Type of `a2`
 * @tparam T_a3 Type of `a3`
 * @tparam T_b1 Type of `b1`
 * @tparam T_b2 Type of `b2`
 * @tparam T_z Type of `z`
 * @param a1 Variable to check
 * @param a2 Variable to check
 * @param a3 Variable to check
 * @param b1 Variable to check
 * @param b2 Variable to check
 * @param z Variable to check
 * @return `true` if `3F2(a1, a2, a3, b1, b2, z)` meets convergence conditions,
 *  and no coefficient is NaN
 */
template <typename T_a1, typename T_a2, typename T_a3, typename T_b1,
          typename T_b2, typename T_z>
inline bool is_3F2_converges(const T_a1& a1, const T_a2& a2, const T_a3& a3,
                             const T_b1& b1, const T_b2& b2, const T_z& z) {
  using std::fabs;
  using std::floor;

  if (!is_not_nan(a1) && !is_not_nan(a2) && !is_not_nan(a3) && !is_not_nan(b1)
      && !is_not_nan(b2) && !is_not_nan(z))
    return false;

  int num_terms = 0;
  bool is_polynomial = false;

  if (is_nonpositive_integer(a1) && fabs(a1) >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs(value_of_rec(a1)));
  }
  if (is_nonpositive_integer(a2) && fabs(a2) >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs(value_of_rec(a2)));
  }
  if (is_nonpositive_integer(a3) && fabs(a3) >= num_terms) {
    is_polynomial = true;
    num_terms = floor(fabs(value_of_rec(a3)));
  }

  if (is_nonpositive_integer(b1)
      && fabs(b1) <= num_terms
      && is_nonpositive_integer(b2)
      && fabs(b2) <= num_terms)
    return false;

  if (!is_polynomial
      && !(fabs(z) < 1.0)
      && !(fabs(z) == 1.0 && b1 + b2 > a1 + a2 + a3))
    return false;

  return true;
}

}  // namespace math
}  // namespace stan
#endif
