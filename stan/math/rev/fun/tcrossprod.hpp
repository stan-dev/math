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

  Eigen::MatrixXd res_val(M.rows(), M.rows());
  arena_matrix<promote_scalar_t<var, T>> res
      = arena_M_val * arena_M_val.transpose();

  reverse_pass_callback([res, arena_M, arena_M_val]() mutable {
    ref_type_t<decltype(res.adj())> adj = res.adj();
    arena_M.adj() += (adj.transpose() + adj) * arena_M_val;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
