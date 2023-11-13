#ifndef STAN_MATH_PRIM_FUN_CHOOSE_HPP
#define STAN_MATH_PRIM_FUN_CHOOSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the binomial coefficient for the specified integer
 * arguments.
 *
 * The binomial coefficient, \f${n \choose k}\f$, read "n choose k", is
 * defined for \f$0 \leq k \leq n\f$ (otherwise return 0) by
 *
 * \f${n \choose k} = \frac{n!}{k! (n-k)!}\f$.
 *
 * @param n total number of objects
 * @param k number of objects chosen
 * @return n choose k or 0 iff k > n
 * @throw std::domain_error if either argument is negative or the
 * result will not fit in an int type
 */
inline int choose(int n, int k) {
  check_nonnegative("choose", "n", n);
  check_nonnegative("choose", "k", k);
  if (k > n) {
    return 0;
  }
  const double choices = boost::math::binomial_coefficient<double>(n, k);
  check_less_or_equal("choose", "n choose k", choices,
                      std::numeric_limits<int>::max());
  return static_cast<int>(std::round(choices));
}

/**
 * Enables the vectorized application of the binomial coefficient function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Binomial coefficient function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto choose(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return choose(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
