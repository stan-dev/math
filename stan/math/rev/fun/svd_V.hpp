#ifndef STAN_MATH_REV_FUN_SVD_V_HPP
#define STAN_MATH_REV_FUN_SVD_V_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Given input matrix m, return matrix V where `m = UDV^{T}`
 *
 * Adjoint update equation comes from Equation (4) in Differentiable Programming
 * Tensor Networks(H. Liao, J. Liu, et al., arXiv:1903.09650).
 *
 * @tparam EigMat type of input matrix
 * @param m MxN input matrix
 * @return Orthogonal matrix V
 */
template <typename EigMat, require_rev_matrix_t<EigMat>* = nullptr>
inline auto svd_V(const EigMat& m) {
  using ret_type = return_var_matrix_t<Eigen::MatrixXd, EigMat>;
  if (unlikely(m.size() == 0)) {
    return ret_type(Eigen::MatrixXd(0, 0));
  }

  const int M = std::min(m.rows(), m.cols());
  auto arena_m = to_arena(m);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      arena_m.val(), Eigen::ComputeThinU | Eigen::ComputeThinV);

  auto arena_D = to_arena(svd.singularValues());

  arena_t<Eigen::MatrixXd> arena_Fm(M, M);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (j == i) {
        arena_Fm(i, j) = 0.0;
      } else {
        arena_Fm(i, j)
            = 1.0 / (arena_D[j] - arena_D[i]) - 1.0 / (arena_D[i] + arena_D[j]);
      }
    }
  }

  auto arena_U = to_arena(svd.matrixU());
  arena_t<ret_type> arena_V = svd.matrixV();

  reverse_pass_callback([arena_m, arena_U, arena_D, arena_V,
                         arena_Fm]() mutable {
    Eigen::MatrixXd VTVadj = arena_V.val_op().transpose() * arena_V.adj_op();
    arena_m.adj()
        += 0.5 * arena_U
               * (arena_Fm.array() * (VTVadj - VTVadj.transpose()).array())
                     .matrix()
               * arena_V.val_op().transpose()
           + arena_U * arena_D.asDiagonal().inverse()
                 * arena_V.adj_op().transpose()
                 * (Eigen::MatrixXd::Identity(arena_m.cols(), arena_m.cols())
                    - arena_V.val_op() * arena_V.val_op().transpose());
  });

  return ret_type(arena_V);
}

}  // namespace math
}  // namespace stan

#endif
