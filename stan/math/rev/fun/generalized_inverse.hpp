#ifndef STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/rev/fun/tcrossprod.hpp>
#include <stan/math/prim/fun/crossprod.hpp>
#include <stan/math/prim/fun/inverse_spd.hpp>
#include <stan/math/rev/fun/mdivide_left_spd.hpp>
#include <stan/math/prim/fun/mdivide_right_spd.hpp>
#include <stan/math/rev/fun/inverse.hpp>

namespace stan {
namespace math {

/*
 * Reverse mode differentiation algorithm reference:
 *
 * The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems
 * Whose Variables Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM
 * Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432
 *
 * Equation 4.12 in the paper
 *
 *  See also
 *  http://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse
 *
 *  Implementation is based on
 *   Title: Uncertainties Python Package
 *    Author: Eric O. LEBIGOT
 *    Date: 2020
 *    Availability:
 * https://github.com/lebigot/uncertainties/blob/master/uncertainties/unumpy/core.py
 */
template <typename VarMat, require_rev_matrix_t<VarMat>* = nullptr>
inline auto generalized_inverse(const VarMat& G) {
  using value_t = value_type_t<VarMat>;
  using ret_type = promote_var_matrix_t<VarMat, VarMat>;

  if (G.size() == 0)
    return ret_type(G);

  if (G.rows() == G.cols())
    return inverse(G);

  if (G.rows() < G.cols()) {
    arena_t<VarMat> G_arena(G);
    auto A_spd = tcrossprod(G_arena.val_op());
    arena_t<VarMat> inv_G(
        transpose(mdivide_left_spd(A_spd.val_op(), G_arena.val_op())));

    auto PG = to_arena(-G_arena.val_op() * inv_G.val_op());
    PG.diagonal().array() += 1.0;

    auto GP = to_arena(-inv_G.val_op() * G_arena.val_op());
    GP.diagonal().array() += 1.0;

    reverse_pass_callback([G_arena, inv_G, GP, PG]() mutable {
      G_arena.adj() += -inv_G.val_op() * G_arena.adj_op() * inv_G.val_op();
      G_arena.adj() += tcrossprod(inv_G.val_op()) * G_arena.adj_op().transpose()
                       * PG.val_op();
      G_arena.adj() += GP.val_op() * G_arena.adj().transpose()
                       * crossprod(inv_G.val_op());
    });
    return ret_type(inv_G);
  } else {
    arena_t<VarMat> G_arena(G);
    auto A_spd = crossprod(G_arena.val_op());
    arena_t<VarMat> inv_G(
        transpose(mdivide_right_spd(G_arena.val_op(), A_spd.val_op())));

    auto PG = to_arena(-G_arena.val_op() * inv_G.val_op());
    PG.diagonal().array() += 1.0;

    auto GP = to_arena(-inv_G.val_op() * G_arena.val_op());
    GP.diagonal().array() += 1.0;

    reverse_pass_callback([G_arena, inv_G, GP, PG]() mutable {
      G_arena.adj() += -inv_G.val_op() * G_arena.adj_op() * inv_G.val_op();
      G_arena.adj() += tcrossprod(inv_G.val_op()) * G_arena.adj_op().transpose()
                       * PG.val_op();
      G_arena.adj() += GP.val_op() * G_arena.adj().transpose()
                       * crossprod(inv_G.val_op());
    });
    return ret_type(inv_G);
  }
}

/*
 * Reverse mode differentiation algorithm reference:
 *
 * The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems
 * Whose Variables Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM
 * Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432
 *
 * Equation 4.12 in the paper
 *  See also
 *  http://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse
 *
 *  Implementation is based on
 *   Title: Uncertainties Python Package
 *    Author: Eric O. LEBIGOT
 *    Date: 2020
 *    Availability:
 * https://github.com/lebigot/uncertainties/blob/master/uncertainties/unumpy/core.py
 */
template <typename VarMat, require_rev_matrix_t<VarMat>* = nullptr>
inline auto generalized_inverse(const VarMat& G, const double a) {
  using value_t = value_type_t<VarMat>;
  using ret_type = promote_var_matrix_t<VarMat, VarMat>;

  if (G.size() == 0)
    return G;

  if (G.rows() == G.cols())
    return inverse(G);

  if (G.rows() < G.cols()) {
    arena_t<VarMat> G_arena(G);
    auto A_spd = tcrossprod(G_arena.val_op());
    A_spd.diagonal().array() += a;
    arena_t<VarMat> inv_G(transpose(mdivide_left_spd(A_spd, G_arena.val_op())));

    auto PG = to_arena(-G_arena.val_op() * inv_G.val_op());
    PG.diagonal().array() += 1.0;

    auto GP = to_arena(-inv_G.val_op() * G_arena.val_op());
    GP.diagonal().array() += 1.0;

    reverse_pass_callback([G_arena, inv_G, GP, PG]() mutable {
      G_arena.adj() += -inv_G.val_op() * G_arena.adj_op() * inv_G.val_op();
      G_arena.adj() += tcrossprod(inv_G.val_op()) * G_arena.adj_op().transpose()
                       * PG.val_op();
      G_arena.adj() += GP.val_op() * G_arena.adj().transpose()
                       * crossprod(inv_G.val_op());
    });
    return ret_type(inv_G);
  } else {
    arena_t<VarMat> G_arena(G);
    auto A_spd = crossprod(G_arena.val_op());
    A_spd.diagonal().array() += a;
    arena_t<VarMat> inv_G(
        transpose(mdivide_right_spd(G_arena.val_op(), A_spd)));

    auto PG = to_arena(-G_arena.val_op() * inv_G.val_op());
    PG.diagonal().array() += 1.0;

    auto GP = to_arena(-inv_G.val_op() * G_arena.val_op());
    GP.diagonal().array() += 1.0;

    reverse_pass_callback([G_arena, inv_G, GP, PG]() mutable {
      G_arena.adj() += -inv_G.val_op() * G_arena.adj_op() * inv_G.val_op();
      G_arena.adj() += tcrossprod(inv_G.val_op()) * G_arena.adj_op().transpose()
                       * PG.val_op();
      G_arena.adj() += GP.val_op() * G_arena.adj().transpose()
                       * crossprod(inv_G.val_op());
    });
    return ret_type(inv_G);
  }
}

}  // namespace math
}  // namespace stan

#endif
