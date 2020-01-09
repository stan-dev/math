#ifndef STAN_MATH_PRIM_FUN_REP_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_REP_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T>
inline Eigen::Matrix<return_type_t<T>, Eigen::Dynamic, 1> rep_vector(const T& x,
                                                                     int n) {
  check_nonnegative("rep_vector", "n", n);
  return Eigen::Matrix<return_type_t<T>, Eigen::Dynamic, 1>::Constant(n, x);
}

}  // namespace math
}  // namespace stan

#endif
