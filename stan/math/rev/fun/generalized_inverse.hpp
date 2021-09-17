#ifndef STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/fun/add_diag.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/inverse.hpp>

namespace stan {
namespace math {

namespace internal {
/*
 * Reverse mode specialization of calculating the generalized inverse of a
 * matrix.
 * <ul><li> Golub, G.H. and Pereyra, V. The Differentiation of Pseudo-Inverses
 * and Nonlinear Least Squares Problems Whose Variables Separate. <i>SIAM
 * Journal on Numerical Analysis</i>, Vol. 10, No. 2 (Apr., 1973), pp.
 * 413-432</li></ul>
 */
template <typename T1, typename T2>
inline auto generalized_inverse_lambda(T1& G_arena, T2& inv_G) {
  return [G_arena, inv_G]() mutable {
    G_arena.adj()
        += -(inv_G.val_op().transpose() * inv_G.adj_op()
             * inv_G.val_op().transpose())
           + (-G_arena.val_op() * inv_G.val_op()
              + Eigen::MatrixXd::Identity(G_arena.rows(), inv_G.cols()))
                 * inv_G.adj_op().transpose() * inv_G.val_op()
                 * inv_G.val_op().transpose()
           + inv_G.val_op().transpose() * inv_G.val_op()
                 * inv_G.adj_op().transpose()
                 * (-inv_G.val_op() * G_arena.val_op()
                    + Eigen::MatrixXd::Identity(inv_G.rows(), G_arena.cols()));
  };
}
}  // namespace internal

/*
 * Reverse mode specialization of calculating the generalized inverse of a
 * matrix.
 *
 * @param G specified matrix
 * @return Generalized inverse of the matrix (an empty matrix if the specified
 * matrix has size zero).
 *
 * @note For the derivatives of this function to exist the matrix must be
 * of constant rank.
 * Reverse mode differentiation algorithm reference:
 *
 * <ul><li> Golub, G.H. and Pereyra, V. The Differentiation of Pseudo-Inverses
 * and Nonlinear Least Squares Problems Whose Variables Separate. <i>SIAM
 * Journal on Numerical Analysis</i>, Vol. 10, No. 2 (Apr., 1973), pp.
 * 413-432</li></ul>
 *
 * Equation 4.12 in the paper
 *
 *  See also
 *  http://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse
 *
 */
template <typename VarMat, require_rev_matrix_t<VarMat>* = nullptr>
inline auto generalized_inverse(const VarMat& G) {
  using ret_type = promote_var_matrix_t<VarMat, VarMat>;

  if (G.size() == 0)
    return ret_type(G);

  if (G.rows() == G.cols()) {
    arena_t<VarMat> G_arena(G);
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>
        complete_ortho_decomp_G
        = G_arena.val().completeOrthogonalDecomposition();
    if (!(complete_ortho_decomp_G.rank() < G.rows())) {
      return ret_type(inverse(G_arena));
    } else {
      arena_t<ret_type> inv_G(complete_ortho_decomp_G.pseudoInverse());
      reverse_pass_callback(
          internal::generalized_inverse_lambda(G_arena, inv_G));
      return ret_type(inv_G);
    }
  } else if (G.rows() < G.cols()) {
    arena_t<VarMat> G_arena(G);
    arena_t<ret_type> inv_G((G_arena.val_op() * G_arena.val_op().transpose())
                                .ldlt()
                                .solve(G_arena.val_op())
                                .transpose());
    reverse_pass_callback(internal::generalized_inverse_lambda(G_arena, inv_G));
    return ret_type(inv_G);
  } else {
    arena_t<VarMat> G_arena(G);
    arena_t<ret_type> inv_G((G_arena.val_op().transpose() * G_arena.val_op())
                                .ldlt()
                                .solve(G_arena.val_op().transpose()));
    reverse_pass_callback(internal::generalized_inverse_lambda(G_arena, inv_G));
    return ret_type(inv_G);
  }
}

}  // namespace math
}  // namespace stan
#endif
