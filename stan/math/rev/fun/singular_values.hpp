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
 * Return the singular values of the specified symmetric matrix.
 *
 * For MxN input matrix, M must be <= N
 *
 * @tparam EigMat type of input matrix
 * @param m MxN input matrix, M <= N
 * @return Singular values of matrix
 */
template <typename EigMat,
	  require_eigen_vt<EigMat>* = nullptr,
	  require_not_vt_var<EigMat>* = nullptr>
inline auto singular_values(const EigMat& m) {
  using ret_type = promote_scalar_t<var, Eigen::VectorXd>;
  check_nonzero_size("singular_values", "m", m);

  const int M = m.rows();
  auto arena_m = to_arena(m);
  
  Eigen::MatrixXd m_val = value_of(arena_m);
  Eigen::JacobiSVD<EigMat> svd(m_val, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> D = Eigen::MatrixXd::Zero(M,M);
  D.diagonal() = svd.singularValues();
  arena_t<ret_type> sigular_values = svd.singularValues();

  auto U = svd.matrixU();
  auto V = svd.matrixU();

  // equation (4) from "Differentiable Programming Tensor Networks"
  Eigen::MatrixXd Fp(M,M);
  Eigen::MatrixXd Fm(M,M);

  for(int i = 0; i < M; i++) {
    for(int j = 0; i < M; j++) {
      if(j == i) {
	Fp[i, j] = 0.0;
	Fm[i, j] = 0.0;
      } else {
	Fp[i, j] = 1.0 / (D[j,j] - D[i,i]) + 1.0 / (D[i,i] + D[j,j]);
	Fm[i, j] = 1.0 / (D[j,j] - D[i,i]) - 1.0 / (D[i,i] + D[j,j]);
      }
    }
  }

  reverse_pass_callback([arena_m, U, V]() mutable {
    arena_m.adj() += 0.5 * U * (Fp * (U.transpose() * U.adj() - t(U.adj()) * U)) * V.transpose() // adjU contributions
      +(Eigen::MatrixXd::Identity(M, M) - U * U.transpose()) * U.adj() * D.inverse() * V.transpose()
      + U * D.adj() * V.transpose() // adjD contributions
      + 0.5 * U * (Fm * (V.transpose() * V.adj() - V.adj().transpose() * V)) * V.transpose() // adjV contributions
      + U * D.inverse() * V.adj().transpose() * (Eigen::MatrixXd::Identity(M, M) - V * V.transpose());
  });
  return ret_type(sigularvalues);
}
}  // namespace math
}  // namespace stan

#endif
