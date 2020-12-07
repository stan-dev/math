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

  if (G.size() == 0) 
    return G;
  
  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols()) 
    return inverse(G);
  
  if (n < m) {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = tcrossprod(G_arena.val());
    arena_t<EigMat> inv_G(transpose(mdivide_left_spd(A_spd, G_arena.val())));
    auto aP = (1 - G_arena.val() * inv_G.transpose());
    auto Pa = (1 - inv_G * G_arena.val());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      G_arena.adj() += -inv_G * G_arena.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G_arena.adj().transpose() * aP;
      G_arena.adj() += Pa * G_arena.adj().transpose() * crossprod(inv_G());
    });
    return inv_G;
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = crossprod(G_arena.val());
    arena_t<EigMat> inv_G(transpose(mdivide_right_spd(G_arena.val(), A_spd)));
    auto aP = (1 - G_arena.val() * inv_G.transpose());
    auto Pa = (1 - inv_G * G_arena.val());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      G_arena.adj() += -inv_G * G_arena.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G_arena.adj().transpose() * aP;
      G_arena.adj() += Pa * G_arena.adj().transpose() * crossprod(inv_G());
    });
    return inv_G;
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
  
  if (G.size() == 0) 
    return G;

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols()) 
    return inverse(G);

  if (n < m) {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = tcrossprod(G_arena.val());
    A_spd.diagonal().array() += a;

    arena_t<EigMat> inv_G(transpose(mdivide_left_spd(A_spd, G_arena.val())));
    auto aP = (1 - G_arena.val() * inv_G.transpose());
    auto Pa = (1 - inv_G * G_arena.val());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      G_arena.adj() += -inv_G * G_arena.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G_arena.adj().transpose() * aP;
      G_arena.adj() += Pa * G_arena.adj().transpose() * crossprod(inv_G);
    });
    return inv_G;
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    auto A_spd = crossprod(G_arena.val());
    A_spd.diagonal().array() += a;

    arena_t<EigMat> inv_G(transpose(mdivide_right_spd(G_arena.val(), A_spd)));
    auto aP = (1 - G * inv_G.transpose());
    auto Pa = (1 - inv_G * G_arena.val());
    reverse_pass_callback([G_arena, inv_G, aP, Pa]() mutable {
      G_arena.adj() += -inv_G * G_arena.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G_arena.adj().transpose() * aP;
      G_arena.adj() += Pa * G_arena.adj().transpose() * crossprod(inv_G());
    });
    return inv_G;
  }
}

}  // namespace math
}  // namespace stan

#endif
