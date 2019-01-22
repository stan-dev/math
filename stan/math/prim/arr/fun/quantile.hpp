#ifndef STAN_MATH_PRIM_ARR_FUN_QUANTILE_HPP
#define STAN_MATH_PRIM_ARR_FUN_QUANTILE_HPP

#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/arr/fun/sort_asc.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the the q-th quantile of the specified sequence of values.
 *
 * @param[in] x array of specified values
 * @param[in] p quantile to calculate, which must be between 0 and 1
 * @return The p-th quantile of the vector values
 */
inline double quantile(const std::vector<double>& x, double p) {
  using std::ceil;
  using std::floor;

  check_greater_or_equal("quantile", "p", p, 0.0);
  check_less_or_equal("quantile", "p", p, 1.0);

  std::vector<double> sorted = sort_asc(x);

  double index = p * (sorted.size() - 1);
  size_t ceiled = static_cast<size_t>(ceil(index));
  size_t floored = static_cast<size_t>(floor(index));

  if (ceiled == floored) {
    return sorted[ceiled];
  }

  return (ceiled - index) * sorted[floored]
         + (index - floored) * sorted[ceiled];
}

}  // namespace math
}  // namespace stan

#endif
