#ifndef STAN_MATH_PRIM_FUN_REP_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_REP_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T_ret, typename T,
          require_eigen_col_vector_t<T_ret>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
inline auto rep_vector(const T& x, int n) {
  check_nonnegative("rep_vector", "n", n);
  return T_ret::Constant(n, x);
}
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto rep_vector(const T& x, int n) {
  return rep_vector<Eigen::Matrix<return_type_t<T>, Eigen::Dynamic, 1>>(x, n);
}

}  // namespace math
}  // namespace stan

#endif
