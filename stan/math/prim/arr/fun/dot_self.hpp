#ifndef STAN_MATH_PRIM_ARR_FUN_DOT_SELF_HPP
#define STAN_MATH_PRIM_ARR_FUN_DOT_SELF_HPP

#include <stan/math/prim/meta.hpp>
#include <vector>
#include <numeric>
#include <cstddef>

namespace stan {
namespace math {

inline double dot_self(const std::vector<double>& x) {
  return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
}

}  // namespace math
}  // namespace stan
#endif
