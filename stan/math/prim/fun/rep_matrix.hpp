#ifndef STAN_MATH_PRIM_FUN_REP_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_REP_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto rep_matrix(const T& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
#ifdef USE_STANC3
  return Eigen::Matrix<return_type_t<T>, Eigen::Dynamic,
                       Eigen::Dynamic>::Constant(m, n, x);
#else
  return Eigen::Matrix<return_type_t<T>, Eigen::Dynamic,
                       Eigen::Dynamic>::Constant(m, n, x).eval();
#endif
}

template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
inline auto rep_matrix(const ColVec& v, int n) {
  check_nonnegative("rep_matrix", "rows", n);
#ifdef USE_STANC3
  return v.replicate(1, n);
#else
  return v.replicate(1, n).eval();
#endif
}

template <typename RowVec, require_eigen_row_vector_t<RowVec>* = nullptr>
inline auto rep_matrix(const RowVec& rv, int m) {
  check_nonnegative("rep_matrix", "cols", m);
#ifdef USE_STANC3
  return rv.replicate(m, 1);
#else
  return rv.replicate(m, 1).eval();
#endif
}
}  // namespace math
}  // namespace stan

#endif
