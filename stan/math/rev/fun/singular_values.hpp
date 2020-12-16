#ifndef STAN_MATH_REV_FUN_SINGULAR_VALUES_HPP
#define STAN_MATH_REV_FUN_SINGULAR_VALUES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <Eigen/SVD>

namespace stan {
namespace math {

/**
 * Return the singular values of the specified matrix.
 *
 * Equation (4) from Differentiable Programming Tensor Networks(H. Liao, J. Liu, et al., arXiv:1903.09650)
 * is used for MxN input matrix's adjoint calculation.
 *
 * @tparam EigMat type of input matrix
 * @param m MxN input matrix
 * @return Singular values of matrix
 */
template <typename EigMat,
	  require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
	  require_vt_var<EigMat>* = nullptr>
inline auto singular_values(const EigMat& m) {
  using ret_type = promote_scalar_t<var, Eigen::VectorXd>;
  check_nonzero_size("singular_values", "m", m);

  const int M = to_arena(m.rows());
  auto arena_m = to_arena(m);
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(arena_m.val(),
					Eigen::ComputeFullU | Eigen::ComputeFullV);

  arena_t<ret_type> singular_values = svd.singularValues();

  auto U = to_arena(svd.matrixU());
  auto V = to_arena(svd.matrixV());

  reverse_pass_callback(
      [arena_m, U, singular_values, V]() mutable {
	arena_m.adj() += U.block(0, 0, U.rows(), singular_values.size()) *
	  singular_values.adj().asDiagonal() *
	  V.block(0, 0, V.rows(), singular_values.size()).transpose();
      });

  return ret_type(singular_values);
}

}  // namespace math
}  // namespace stan

#endif
