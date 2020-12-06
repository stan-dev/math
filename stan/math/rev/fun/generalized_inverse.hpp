#ifndef STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_REV_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/cholesky_decompose.hpp>
#include <stan/math/rev/fun/transpose.hpp>
#include <stan/math/rev/fun/tcrossprod.hpp>
#include <stan/math/rev/fun/crossprod.hpp>
#include <stan/math/rev/fun/inverse_spd.hpp>
#include <stan/math/rev/fun/mdivide_left_spd.hpp>
#include <stan/math/rev/fun/mdivide_right_spd.hpp>
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
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline auto generalized_inverse(const EigMat& G) {
  using value_t = value_type_t<EigMat>;
  if (G.size() == 0) {
    return {};
  }

  const auto n = G.rows();
  const auto m = G.cols();

  if (A.rows() == G.cols()) {
    return inverse(G);
  }
  if (n < m) {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    arena_t<EigMat> A_spd = tcrossprod(G_arena.val());
    arena_t<EigMat> inv_G(transpose(mdivide_left_spd(A_spd, G_arena.val())));
    arena_t<EigMat> aP = (1 - G * inv_G.transpose());
    arena_t<EigMat> Pa = (1 - inv_G * G);
    reverse_pass_callback([G_arena, inv_G]() mutable {
      G_arena.adj() += -inv_G * G.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G.adj().transpose() * aP;
      G_arena.adj() += Pa * G.adj().transpose() * crossprod(inv_G());
    });
    return ret;
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    arena_t<EigMat> A_spd = crossprod(G_arena.val());
    arena_t<EigMat> inv_G(transpose(mdivide_right_spd(G_arena.val(), A_spd)));
    arena_t<EigMat> aP = (1 - G * inv_G.transpose());
    arena_t<EigMat> Pa = (1 - inv_G * G);
    reverse_pass_callback([G_arena, inv_G]() mutable {
      G_arena.adj() += -inv_G * G.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G.adj().transpose() * aP;
      G_arena.adj() += Pa * G.adj().transpose() * crossprod(inv_G());
    });
    return ret;
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
template <typename EigMat, typename Scal, require_eigen_t<EigMat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
generalized_inverse(const EigMat& G, const Scal& a) {
  using value_t = value_type_t<EigMat>;
  if (G.size() == 0) {
    return {};
  }

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols()) {
    return inverse(G);
  }

  if (n < m) {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    arena_t<EigMat> A_spd = tcrossprod(G_arena.val());
    A_spd.diagonal().array() += a;
<<<<<<< HEAD
    
    arena_t<EigMat> inv_G(transpose(mdivide_left_spd(A_spd, G_arena.val())));
    arena_t<EigMat> aP = (1 - G * inv_G.transpose());
    arena_t<EigMat> Pa = (1 - inv_G * G);
     reverse_pass_callback([G_arena, inv_G]() mutable {
      G_arena.adj() += -inv_G * G.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G.adj().transpose() * aP;
      G_arena.adj() += Pa * G.adj().transpose() * crossprod(inv_G());
=======

    arena_t<EigMat> inv_A(transpose(mdivide_left_spd(A_spd, A_arena.val())));
    arena_t<EigMat> aP = (1 - A * inv_A.transpose());
    arena_t<EigMat> Pa = (1 - inv_A * A);
    reverse_pass_callback([A_arena, inv_A]() mutable {
      A_arena.adj() += -inv_A * A.adj() * inv_A;
      A_arena.adj() += tcrossprod(inv_A) * A.adj().transpose() * aP;
      A_arena.adj() += Pa * A.adj().transpose() * crossprod(inv_A());
>>>>>>> 9e6e692cb1f69058e1bf293800d9fa7660f3bc7f
    });
    return ret;
  } else {
    arena_t<plain_type_t<EigMat>> G_arena(G);
    arena_t<EigMat> A_spd = crossprod(G_arena.val());
    A_spd.diagonal().array() += a;

    arena_t<EigMat> inv_G(transpose(mdivide_right_spd(G_arena.val(), A_spd));
    arena_t<EigMat> aP = (1 - G * inv_G.transpose());
    arena_t<EigMat> Pa = (1 - inv_G * G);
   reverse_pass_callback([G_arena, inv_G]() mutable {
      G_arena.adj() += -inv_G * G.adj() * inv_G;
      G_arena.adj() += tcrossprod(inv_G) * G.adj().transpose() * aP;
      G_arena.adj() += Pa * G.adj().transpose() * crossprod(inv_G());
    });
    return ret;
  }
}

}  // namespace math
}  // namespace stan

#endif
