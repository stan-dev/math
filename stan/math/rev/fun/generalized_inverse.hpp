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
inline auto generalized_inverse(const EigMat& A) {
  using value_t = value_type_t<EigMat>;
  if (A.size() == 0) {
    return {};
  }

  const auto n = A.rows();
  const auto m = A.cols();

  if (A.rows() == A.cols()) {
    return inverse(A);
  }
  if (n < m) {
    arena_t<plain_type_t<EigMat>> A_arena(A);
    arena_t<EigMat> A_spd = tcrossprod(A_arena.val());
    arena_t<EigMat> inv_A(transpose(mdivide_left_spd(A_spd, A_arena.val())));
    arena_t<EigMat> aP = (1 - A * inv_A.transpose());
    arena_t<EigMat> Pa = (1 - inv_A * A);
    reverse_pass_callback([A_arena, inv_A]() mutable {
      A_arena.adj() += -inv_A * A.adj() * inv_A;
      A_arena.adj() += tcrossprod(inv_A) * A.adj().transpose() * aP;
      A_arena.adj() += Pa * A.adj().transpose() * crossprod(inv_A());
    });
    return ret;
  } else {
    arena_t<plain_type_t<EigMat>> A_arena(A);
    arena_t<EigMat> A_spd = crossprod(A_arena.val());
    arena_t<EigMat> inv_A(transpose(mdivide_right_spd(A_arena.val(), A_spd)));
    arena_t<EigMat> aP = (1 - A * inv_A.transpose());
    arena_t<EigMat> Pa = (1 - inv_A * A);
    reverse_pass_callback([A_arena, inv_A]() mutable {
      A_arena.adj() += -inv_A * A.adj() * inv_A;
      A_arena.adj() += tcrossprod(inv_A) * A.adj().transpose() * aP;
      A_arena.adj() += Pa * A.adj().transpose() * crossprod(inv_A());
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
generalized_inverse(const EigMat& A, const Scal& a) {
  using value_t = value_type_t<EigMat>;
  if (A.size() == 0) {
    return {};
  }

  const auto n = A.rows();
  const auto m = A.cols();

  if (A.rows() == A.cols()) {
    return inverse(A);
  }

  if (n < m) {
    arena_t<plain_type_t<EigMat>> A_arena(A);
    arena_t<EigMat> A_spd = tcrossprod(A_arena.val());
    A_spd.diagonal().array() += a;
    
    arena_t<EigMat> inv_A(transpose(mdivide_left_spd(A_spd, A_arena.val())));
    arena_t<EigMat> aP = (1 - A * inv_A.transpose());
    arena_t<EigMat> Pa = (1 - inv_A * A);
    reverse_pass_callback([A_arena, inv_A]() mutable {
      A_arena.adj() += -inv_A * A.adj() * inv_A;
      A_arena.adj() += tcrossprod(inv_A) * A.adj().transpose() * aP;
      A_arena.adj() += Pa * A.adj().transpose() * crossprod(inv_A());
    });
    return ret;
  } else {
    arena_t<plain_type_t<EigMat>> A_arena(A);
    arena_t<EigMat> A_spd = crossprod(A_arena.val());
    A_spd.diagonal().array() += a;

    arena_t<plain_type_t<EigMat>> A_arena(A);
    arena_t<EigMat> inv_A(transpose(mdivide_right_spd(A_arena.val(), A_spd));
    arena_t<EigMat> aP = (1 - A * inv_A.transpose());
    arena_t<EigMat> Pa = (1 - inv_A * A);
    reverse_pass_callback([A_arena, inv_A]() mutable {
      A_arena.adj() += -inv_A * A.adj() * inv_A;
      A_arena.adj() += tcrossprod(inv_A) * A.adj().transpose() * aP;
      A_arena.adj() += Pa * A.adj().transpose() * crossprod(inv_A());
    });
    return ret;
  }
}

}  // namespace math
}  // namespace stan

#endif
