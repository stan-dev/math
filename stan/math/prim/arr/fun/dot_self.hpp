#ifndef STAN_MATH_PRIM_ARR_FUN_DOT_SELF_HPP
#define STAN_MATH_PRIM_ARR_FUN_DOT_SELF_HPP

#include <cstddef>
#include <vector>

namespace stan {
namespace math {

inline double dot_self(const std::vector<double>& x) {
  double sum = 0.0;
  for (size_t i = 0; i < x.size(); ++i) sum += x[i] * x[i];
  return sum;
}

}  // namespace math
}  // namespace stan
#endif
