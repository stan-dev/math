#ifndef STAN_MATH_REV_FUN_SVD_HPP
#define STAN_MATH_REV_FUN_SVD_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/svd.hpp>

namespace stan {
namespace math {

/**
 * Given input matrix m, return the singular value decomposition (U,D,V)
 * such that `m = U*diag(D)*V^{T}`
 *
 * Adjoint update equation comes from Equation (4) in Differentiable Programming
 * Tensor Networks(H. Liao, J. Liu, et al., arXiv:1903.09650).
 *
 * @tparam EigMat type of input matrix
 * @param m MxN input matrix
 * @return a tuple (U,D,V) where U is an orthogonal matrix, D a vector of
 * singular values (in decreasing order), and V an orthogonal matrix
 */
template <typename EigMat, require_rev_matrix_t<EigMat>* = nullptr>
inline auto svd(const EigMat& m) {
  using mat_ret_type = return_var_matrix_t<Eigen::MatrixXd, EigMat>;
  using vec_ret_type = return_var_matrix_t<Eigen::VectorXd, EigMat>;

  if (unlikely(m.size() == 0)) {
    return std::make_tuple(mat_ret_type(Eigen::MatrixXd(0, 0)),
                           vec_ret_type(Eigen::VectorXd(0, 1)),
                           mat_ret_type(Eigen::MatrixXd(0, 0)));
  }

  const int M = std::min(m.rows(), m.cols());
  auto arena_m = to_arena(m);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      arena_m.val(), Eigen::ComputeThinU | Eigen::ComputeThinV);

  auto singular_values_d = svd.singularValues();

  arena_t<Eigen::MatrixXd> arena_Fp(M, M);
  arena_t<Eigen::MatrixXd> arena_Fm(M, M);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      double a = 1.0 / (singular_values_d[j] - singular_values_d[i]);
      double b = 1.0 / (singular_values_d[i] + singular_values_d[j]);
      arena_Fp(i, j) = a + b;
      arena_Fm(i, j) = a - b;
    }
  }
  arena_Fp.diagonal().setZero();
  arena_Fm.diagonal().setZero();

  arena_t<vec_ret_type> singular_values = singular_values_d;
  arena_t<mat_ret_type> arena_U = svd.matrixU();
  arena_t<mat_ret_type> arena_V = svd.matrixV();

  reverse_pass_callback([arena_m, arena_U, singular_values, arena_V, arena_Fp,
                         arena_Fm]() mutable {
    // SVD-U reverse mode
    Eigen::MatrixXd UUadjT = arena_U.val_op().transpose() * arena_U.adj_op();
    auto u_adj
        = .5 * arena_U.val_op()
              * (arena_Fp.array() * (UUadjT - UUadjT.transpose()).array())
                    .matrix()
              * arena_V.val_op().transpose()
          + (Eigen::MatrixXd::Identity(arena_m.rows(), arena_m.rows())
             - arena_U.val_op() * arena_U.val_op().transpose())
                * arena_U.adj_op()
                * singular_values.val_op().asDiagonal().inverse()
                * arena_V.val_op().transpose();
    // Singular values reverse mode
    auto d_adj = arena_U.val_op() * singular_values.adj().asDiagonal()
                 * arena_V.val_op().transpose();
    // SVD-V reverse mode
    Eigen::MatrixXd VTVadj = arena_V.val_op().transpose() * arena_V.adj_op();
    auto v_adj
        = 0.5 * arena_U.val_op()
              * (arena_Fm.array() * (VTVadj - VTVadj.transpose()).array())
                    .matrix()
              * arena_V.val_op().transpose()
          + arena_U.val_op() * singular_values.val_op().asDiagonal().inverse()
                * arena_V.adj_op().transpose()
                * (Eigen::MatrixXd::Identity(arena_m.cols(), arena_m.cols())
                   - arena_V.val_op() * arena_V.val_op().transpose());

    arena_m.adj() += u_adj + d_adj + v_adj;
  });

  return std::make_tuple(mat_ret_type(arena_U), vec_ret_type(singular_values),
                         mat_ret_type(arena_V));
}

}  // namespace math
}  // namespace stan

#endif
