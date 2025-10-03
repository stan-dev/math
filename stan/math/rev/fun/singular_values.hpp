#ifndef STAN_MATH_REV_FUN_SINGULAR_VALUES_HPP
#define STAN_MATH_REV_FUN_SINGULAR_VALUES_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>

namespace stan {
namespace math {

/**
 * Return the singular values of the specified matrix.
 *
 * Adjoint update equation comes from Equation (4) in Differentiable Programming
 * Tensor Networks(H. Liao, J. Liu, et al., arXiv:1903.09650).
 *
 * @tparam EigMat type of input matrix
 * @param m MxN input matrix
 * @return Singular values of matrix
 */
template <typename EigMat, require_rev_matrix_t<EigMat>* = nullptr>
inline auto singular_values(const EigMat& m) {
  using ret_type = return_var_matrix_t<Eigen::VectorXd, EigMat>;
  if (unlikely(m.size() == 0)) {
    return ret_type(Eigen::VectorXd(0));
  }

  auto arena_m = to_arena(m);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      arena_m.val(), Eigen::ComputeThinU | Eigen::ComputeThinV);

  arena_t<ret_type> singular_values = svd.singularValues();

  auto arena_U = to_arena(svd.matrixU());
  auto arena_V = to_arena(svd.matrixV());

  reverse_pass_callback([arena_m, arena_U, singular_values, arena_V]() mutable {
    arena_m.adj()
        += arena_U * singular_values.adj().asDiagonal() * arena_V.transpose();
  });

  return ret_type(singular_values);
}

}  // namespace math
}  // namespace stan

#endif
