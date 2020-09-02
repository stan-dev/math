#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline Eigen::Matrix<var, T::RowsAtCompileTime, T::RowsAtCompileTime>
tcrossprod(const T& M) {
  if (M.rows() == 0) {
    return {};
  }

  arena_matrix<promote_scalar_t<var, T>> arena_M = M;
  arena_matrix<promote_scalar_t<double, T>> arena_M_val = value_of(arena_M);

  Eigen::Matrix<double, T::RowsAtCompileTime, T::RowsAtCompileTime> res_val(M.rows(), M.rows());
  res_val.setZero().template selfadjointView<Eigen::Upper>().rankUpdate(arena_M_val);

  arena_matrix<Eigen::Matrix<var, T::RowsAtCompileTime, T::RowsAtCompileTime>> res(M.rows(), M.rows());

  for(size_t j = 0; j < res.cols(); ++j) {
    for(size_t i = 0; i < j; ++i) {
      res.coeffRef(i, j) = res.coeffRef(j, i) = res_val.coeff(i, j);
    }
    res.coeffRef(j, j) = res_val.coeff(j, j);
  }

  reverse_pass_callback([res, arena_M, arena_M_val]() mutable {
    Eigen::MatrixXd adj = res.adj();
    for(size_t i = 0; i < adj.cols(); ++i)
      adj(i, i) *= 2.0;
    arena_M.adj() += adj * arena_M_val;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
