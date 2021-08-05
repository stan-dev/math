#ifndef STAN_MATH_PRIM_FUN_REP_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_REP_ROW_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T>
inline auto rep_row_vector(const T& x, int m) {
  check_nonnegative("rep_row_vector", "m", m);
  return Eigen::Matrix<return_type_t<T>, 1, Eigen::Dynamic>::Constant(m, x);
}

}  // namespace math
}  // namespace stan

#endif
