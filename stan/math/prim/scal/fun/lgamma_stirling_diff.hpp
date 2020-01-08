#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP

#include <cmath>

namespace stan {
namespace math {

double lgamma_stirling_diff(const double x) { return std::log1p(1 / (12 * x)); }

}  // namespace math
}  // namespace stan

#endif
