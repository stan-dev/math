#ifndef STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/fun/add_diag.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/inverse.hpp>
#include <stan/math/prim/fun/crossprod.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>

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
    auto ones1
        = Eigen::VectorXd::Constant(std::min(G_arena.rows(), inv_G.cols()), 1);
    auto ones2
        = Eigen::VectorXd::Constant(std::min(inv_G.rows(), G_arena.cols()), 1);
    G_arena.adj()
        += -(inv_G.val_op().transpose() * inv_G.adj_op()
             * inv_G.val_op().transpose())
           + add_diag(-G_arena.val_op() * inv_G.val_op(), ones1)
                 * inv_G.adj_op().transpose() * inv_G.val_op()
                 * inv_G.val_op().transpose()
           + inv_G.val_op().transpose() * inv_G.val_op()
                 * inv_G.adj_op().transpose()
                 * add_diag(-inv_G.val_op() * G_arena.val_op(), ones2);
  };
}

template <typename AMat>
inline auto cholesky_low_rank_decomposition(const AMat& A) {
    // if (A.size() == 0)
    // return {};

    int n = A.rows();
    auto dA = A.diagonal();
    double tol = 1e-12;
    plain_type_t<AMat> L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
    int r = 0;

    for (int i = 0; i < n; i++) {
      if (r == 0) {
        L.block(i, r, n - i, 1) = A.block(i, i, n - i, 1);
      } else {
        L.block(i, r , n - i, 1) = A.block(i, i, n - i, 1) - (L.block(i, 0, n - i, r) * L.block(i, 0, 1, r).transpose());
      }

      if (L(i, r) > tol ) {
        L(i, r) = sqrt(L(i, r)) ;
        if (i < n - 1) {
             L.block(i + 1, r, n - (i + 1), 1) = (L.block(i + 1, r, n - (i + 1), 1).array() / L(i, r)).matrix();
        }
      }  else {
        r -= 1;
      }
      r += 1;
    }
  return L.block(0, 0, n, r);
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

  const int n = std::min(G.rows(), G.cols());
  const bool transpose_bool = G.rows() == n ? true : false;
  auto A = transpose_bool ? tcrossprod(G) : crossprod(G);

// note: L may not be square. So L * M ops don't cancel.
  auto arena_L = to_arena(internal::cholesky_low_rank_decomposition(A));
  arena_t<VarMat> arena_M = inverse(crossprod(arena_L.val_op()));

  if (transpose_bool) {
    arena_t<VarMat> G_arena(G);
    arena_t<ret_type> inv_G(G_arena.val_op().transpose() * arena_L.val_op() * arena_M.val_op() * arena_M.val_op() * arena_L.val_op().transpose());
    reverse_pass_callback(internal::generalized_inverse_lambda(G_arena, inv_G));
    return ret_type(inv_G);
  } else {
    arena_t<VarMat> G_arena(G);
     arena_t<ret_type> inv_G(arena_L.val_op() * arena_M.val_op() * arena_M.val_op() * arena_L.val_op().transpose() * G_arena.val_op().transpose());
    reverse_pass_callback(internal::generalized_inverse_lambda(G_arena, inv_G));
    return ret_type(inv_G);
  }

}

}  // namespace math
}  // namespace stan
#endif
