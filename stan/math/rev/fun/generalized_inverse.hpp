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

/**
 * Reverse mode differentiation algorithm reference:
 *
 * The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems
 * Whose Variables Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM
 * Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432
 *
 * Equation 4.12 in the paper
 */
template <typename EigMat, require_rev_matrix_t<EigMat>* = nullptr>
inline auto generalized_inverse(const EigMat& G) {
  using value_t = value_type_t<EigMat>;
  using ret_type = promote_var_matrix_t<EigMat, EigMat>;
  if (G.size() == 0)
    return ret_type(G);

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols())
    return inverse(G);

  if (G.rows() < G.cols()) {
    arena_t<EigMat> G_arena(G);
    arena_t<EigMat> inv_G((G_arena.val_op() * G_arena.val_op().transpose()).llt().solve(G_arena.val_op()));
    reverse_pass_callback([G_arena, inv_G]() mutable {
//      (2, 3) += [(3, 2) * (2, 3) * (3, 2)](3, 2)
      G_arena.adj() += (-inv_G.val_op() * G_arena.adj_op() * inv_G.val_op());
//      (2, 3) += [(3, 2) * (2, 3) * (3, 2) * (2, 2)](3, 2)
      G_arena.adj() += (inv_G.val_op() * inv_G.val_op().transpose() * G_arena.adj_op().transpose() * (1 - (G_arena.val_op() * inv_G.val_op()).array()).matrix());
//      (2, 3) += [(3, 3) * (3, 2) * (2, 3) * (3, 2)
      G_arena.adj() += ((1 - (inv_G.val_op() * G_arena.val_op()).array()).matrix() * G_arena.adj_op().transpose() * inv_G.val_op().transpose() * inv_G.val_op());
    });
    return ret_type(inv_G);
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = crossprod(G_arena.val_op());
    arena_t<EigMat> inv_G(mdivide_right_spd(G_arena.val_op(), A_spd).transpose());
    auto aP = to_arena((1 - (G_arena.val_op() * inv_G.val_op()).array()).matrix());
    auto Pa = to_arena((1 - (inv_G.val_op() * G_arena.val_op()).array()).matrix());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2)](3, 2)
            G_arena.adj() += (-inv_G.val_op() * G_arena.adj_op() * inv_G.val_op()).transpose();
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2) * (2, 2)](3, 2)
            G_arena.adj() += (inv_G.val_op() * inv_G.val_op().transpose() * G_arena.adj_op().transpose() * aP).transpose();
      //      (2, 3) += [(3, 3) * (3, 2) * (2, 3) * (3, 2)
            G_arena.adj() += (Pa * G_arena.adj_op().transpose() * inv_G.val_op().transpose() * inv_G.val_op()).transpose();
    });
    return ret_type(inv_G);
  }
}

/**
 * Reverse mode differentiation algorithm reference:
 *
 * The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems
 * Whose Variables Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM
 * Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432
 *
 * Equation 4.12 in the paper
 */
template <typename EigMat, require_rev_matrix_t<EigMat>* = nullptr>
inline auto generalized_inverse(const EigMat& G, const double a) {
  using value_t = value_type_t<EigMat>;
  using ret_type = promote_var_matrix_t<EigMat, EigMat>;

  if (G.size() == 0)
    return ret_type(G);

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols())
    return inverse(G);

  if (n < m) {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = tcrossprod(G_arena.val_op());
    A_spd.diagonal().array() += a;

    arena_t<EigMat> inv_G(mdivide_left_spd(A_spd, G_arena.val_op()).transpose());
    auto aP = to_arena((1 - (G_arena.val_op() * inv_G.val_op()).array()).matrix());
    auto Pa = to_arena((1 - (inv_G.val_op() * G_arena.val_op()).array()).matrix());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2)](3, 2)
            G_arena.adj() += (-inv_G.val_op() * G_arena.adj_op() * inv_G.val_op()).transpose();
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2) * (2, 2)](3, 2)
            G_arena.adj() += (inv_G.val_op() * inv_G.val_op().transpose() * G_arena.adj_op().transpose() * aP).transpose();
      //      (2, 3) += [(3, 3) * (3, 2) * (2, 3) * (3, 2)
            G_arena.adj() += (Pa * G_arena.adj_op().transpose() * inv_G.val_op().transpose() * inv_G.val_op()).transpose();
    });
    return ret_type(inv_G);
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = crossprod(G_arena.val_op());
    A_spd.diagonal().array() += a;

    arena_t<EigMat> inv_G(mdivide_right_spd(G_arena.val_op(), A_spd).transpose());
    auto aP = to_arena((1 - (G_arena.val_op() * inv_G.val_op()).array()).matrix());
    auto Pa = to_arena((1 - (inv_G.val_op() * G_arena.val_op()).array()).matrix());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2)](3, 2)
            G_arena.adj() += (-inv_G.val_op() * G_arena.adj_op() * inv_G.val_op()).transpose();
      //      (2, 3) += [(3, 2) * (2, 3) * (3, 2) * (2, 2)](3, 2)
            G_arena.adj() += (inv_G.val_op() * inv_G.val_op().transpose() * G_arena.adj_op().transpose() * aP).transpose();
      //      (2, 3) += [(3, 3) * (3, 2) * (2, 3) * (3, 2)
            G_arena.adj() += (Pa * G_arena.adj_op().transpose() * inv_G.val_op().transpose() * inv_G.val_op()).transpose();
    });
    return ret_type(inv_G);
  }
}

}  // namespace math
}  // namespace stan

#endif
