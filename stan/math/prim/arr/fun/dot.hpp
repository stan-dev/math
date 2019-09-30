#ifndef STAN_MATH_PRIM_ARR_FUN_DOT_HPP
#define STAN_MATH_PRIM_ARR_FUN_DOT_HPP

#include <stan/math/prim/meta.hpp>
#include <vector>
#include <numeric>
#include <cstddef>

namespace stan {
namespace math {

inline double dot(const std::vector<double>& x, const std::vector<double>& y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

}  // namespace math
}  // namespace stan
#endif
