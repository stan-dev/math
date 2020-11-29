#ifndef MATH_SINGULAR_VALUES_H
#define MATH_SINGULAR_VALUES_H

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <Eigen/SVD>

namespace stan {
namespace math {

/**
 * Return the singular values of the specified matrix.
 *
 * Equation (4) from Differentiable Programming Tensor Networks
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
  
  Eigen::MatrixXd m_val = value_of(arena_m);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_val, Eigen::ComputeFullU | Eigen::ComputeFullV);

  auto singular_values = to_arena(svd.singularValues());
  auto U = to_arena(svd.matrixU());
  auto V = to_arena(svd.matrixV());
  Eigen::MatrixXd Fp(M,M);
  Eigen::MatrixXd Fm(M,M);
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < M; j++) {
      if(j == i) {
        Fp(i, j) = 0.0;
        Fm(i, j) = 0.0;
      } else {
        Fp(i, j)  = 1.0 / (singular_values[j] - singular_values[i]) + 1.0 / (singular_values[i] + singular_values[j]);
        Fm(i, j) = 1.0 / (singular_values[j] - singular_values[i]) - 1.0 / (singular_values[j] + singular_values[i]);
      }
    }
  }
  auto arena_Fp = to_arena(Fp);
  auto arena_Fm = to_arena(Fm);

  reverse_pass_callback(
      [arena_m, M, arena_Fp, arena_Fm, U, singular_values, V]() mutable {
          arena_m.adj() += 0.5 * U * (arena_Fp * (U.transpose() * U.adj() - U.adj().transpose() * U)) * V.transpose() // adjU contributions
                   +(Eigen::MatrixXd::Identity(M, M) - U * U.transpose()) * U.adj() * singular_values.asDiagonal().inverse() * V.transpose()
                   + U * singular_values.adj().asDiagonal() * V.transpose() // adjD contributions
                   + 0.5 * U * (arena_Fm * (V.transpose() * V.adj() - V.adj().transpose() * V)) * V.transpose() // adjV contributions
                   + U * singular_values.asDiagonal().inverse() * V.adj().transpose() * (Eigen::MatrixXd::Identity(M, M) - V * V.transpose());
      });
  return ret_type(singular_values);
}
}  // namespace math
}  // namespace stan

#endif
