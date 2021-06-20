#ifndef STAN_MATH_PRIM_ERR_CHECK_PFQ_CONVERGES_HPP
#define STAN_MATH_PRIM_ERR_CHECK_PFQ_CONVERGES_HPP

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
template <typename T_p, typename T_q, typename T_z>
inline void check_pFq_converges(const T_p& p, const T_q& q, const T_z& z) {
  using std::fabs;
  using std::floor;

  check_not_nan("check_pFq_converges", "p", p);
  check_not_nan("check_pFq_converges", "q", q);
  check_not_nan("check_pFq_converges", "z", z);

  bool condition_1 = (p.size() < (q.size() + 1));
  bool condition_2 = (p.size() == (q.size() + 1)) && (fabs(z) < 1);
  bool condition_3 = (p.size() > (q.size() + 1)) && (z == 0);

  if (condition_1 || condition_2 || condition_3) {
    return;
  }
 
  std::stringstream msg;
  msg << "hypergeometric function pFq does not meet convergence "
      << "conditions with given arguments. "
      << "p: " << p << ", q: " << q << ", "
      << ", z: " << z;
  throw std::domain_error(msg.str());
}

}  // namespace math
}  // namespace stan
#endif
