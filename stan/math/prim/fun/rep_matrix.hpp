#ifndef STAN_MATH_PRIM_FUN_REP_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_REP_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline Eigen::Matrix<return_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
rep_matrix(const T& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  return Eigen::Matrix<return_type_t<T>, Eigen::Dynamic,
                       Eigen::Dynamic>::Constant(m, n, x);
}

template <typename ColVec, require_eigen_col_vector_t<ColVec>* = nullptr>
inline Eigen::Matrix<value_type_t<ColVec>, Eigen::Dynamic, Eigen::Dynamic>
rep_matrix(const ColVec& v, int n) {
  check_nonnegative("rep_matrix", "rows", n);
  return v.replicate(1, n);
}

template <typename RowVec, require_eigen_row_vector_t<RowVec>* = nullptr>
inline Eigen::Matrix<value_type_t<RowVec>, Eigen::Dynamic, Eigen::Dynamic>
rep_matrix(const RowVec& rv, int m) {
  check_nonnegative("rep_matrix", "cols", m);
  return rv.replicate(m, 1);
}
}  // namespace math
}  // namespace stan

#endif
