#ifndef STAN_MATH_REV_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Compute the trace of an inverse quadratic form premultiplied by a
 * square matrix. This computes
 *       trace(D B^T A^-1 B)
 * where D is a square matrix and the LDLT_factor of A is provided.
 *
 * @tparam Td type of the first matrix
 * @tparam Ta type of matrix in the LDLT_factor
 * @tparam Tb type of the second matrix
 *
 * @param D a square matrix
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <typename Td, typename Ta, typename Tb,
          require_not_col_vector_t<Td>* = nullptr,
          require_all_matrix_t<Td, Ta, Tb>* = nullptr,
          require_any_st_var<Td, Ta, Tb>* = nullptr>
inline var trace_gen_inv_quad_form_ldlt(const Td& D, LDLT_factor<Ta>& A,
                                        const Tb& B) {
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A.matrix(), "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  if (D.size() == 0 || A.matrix().size() == 0) {
    return 0;
  }

  if (!is_constant<Ta>::value && !is_constant<Tb>::value
      && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));
    auto BTAsolveB = to_arena(arena_B.val_op().transpose() * AsolveB);

    var res = (arena_D.val() * BTAsolveB).trace();

    reverse_pass_callback(
        [arena_A, BTAsolveB, AsolveB, arena_B, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().transpose()
                           * AsolveB.transpose();
          arena_B.adj() += C_adj * AsolveB
                           * (arena_D.val_op() + arena_D.val_op().transpose());
          arena_D.adj() += C_adj * BTAsolveB;
        });

    return res;
  } else if (!is_constant<Ta>::value && !is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_D * arena_B.val_op().transpose() * AsolveB).trace();

    reverse_pass_callback([arena_A, AsolveB, arena_B, arena_D, res]() mutable {
      double C_adj = res.adj();

      arena_A.adj()
          -= C_adj * AsolveB * arena_D.transpose() * AsolveB.transpose();
      arena_B.adj() += C_adj * AsolveB * (arena_D + arena_D.transpose());
    });

    return res;
  } else if (!is_constant<Ta>::value && is_constant<Tb>::value
             && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(value_of(B_ref)));
    auto BTAsolveB = to_arena(value_of(B_ref).transpose() * AsolveB);

    var res = (arena_D.val() * BTAsolveB).trace();

    reverse_pass_callback(
        [arena_A, BTAsolveB, AsolveB, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().transpose()
                           * AsolveB.transpose();
          arena_D.adj() += C_adj * BTAsolveB;
        });

    return res;
  } else if (!is_constant<Ta>::value && is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(value_of(B_ref)));

    var res = (arena_D * value_of(B_ref).transpose() * AsolveB).trace();

    reverse_pass_callback([arena_A, AsolveB, arena_D, res]() mutable {
      double C_adj = res.adj();

      arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().transpose()
                       * AsolveB.transpose();
    });

    return res;
  } else if (is_constant<Ta>::value && !is_constant<Tb>::value
             && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));
    auto BTAsolveB = to_arena(arena_B.val_op().transpose() * AsolveB);

    var res = (arena_D.val() * BTAsolveB).trace();

    reverse_pass_callback(
        [BTAsolveB, AsolveB, arena_B, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_B.adj() += C_adj * AsolveB
                           * (arena_D.val_op() + arena_D.val_op().transpose());
          arena_D.adj() += C_adj * BTAsolveB;
        });

    return res;
  } else if (is_constant<Ta>::value && !is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_D * arena_B.val_op().transpose() * AsolveB).trace();

    reverse_pass_callback([AsolveB, arena_B, arena_D, res]() mutable {
      arena_B.adj() += res.adj() * AsolveB * (arena_D + arena_D.transpose());
    });

    return res;
  } else if (is_constant<Ta>::value && is_constant<Tb>::value
             && !is_constant<Td>::value) {
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto BTAsolveB = to_arena(value_of(B_ref).transpose()
                              * A.ldlt().solve(value_of(B_ref)));

    var res = (arena_D.val() * BTAsolveB).trace();

    reverse_pass_callback([BTAsolveB, arena_D, res]() mutable {
      arena_D.adj() += res.adj() * BTAsolveB;
    });

    return res;
  }
}

/**
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       `trace(diag(D) B^T A^-1 B)`
 * where D is the diagonal of a diagonal matrix (`diag(D)` is the diagonal
 * matrix itself) and the LDLT_factor of A is provided.
 *
 * @tparam Td type of the diagonal of first matrix
 * @tparam Ta type of matrix in the LDLT_factor
 * @tparam Tb type of the B matrix
 *
 * @param D diagonal of multiplier
 * @param A LDLT_factor
 * @param B inner term in quadratic form
 * @return trace(diag(D) * B^T * A^-1 * B)
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by diag(D).
 */
template <typename Td, typename Ta, typename Tb,
          require_col_vector_t<Td>* = nullptr,
          require_all_matrix_t<Ta, Tb>* = nullptr,
          require_any_st_var<Td, Ta, Tb>* = nullptr>
inline var trace_gen_inv_quad_form_ldlt(const Td& D, const LDLT_factor<Ta>& A,
                                        const Tb& B) {
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A.matrix(), "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  if (D.size() == 0 || A.matrix().size() == 0) {
    return 0;
  }

  if (!is_constant<Ta>::value && !is_constant<Tb>::value
      && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));
    auto BTAsolveB = to_arena(arena_B.val_op().transpose() * AsolveB);

    var res = (arena_D.val().asDiagonal() * BTAsolveB).trace();

    reverse_pass_callback(
        [arena_A, BTAsolveB, AsolveB, arena_B, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().asDiagonal()
                           * AsolveB.transpose();
          arena_B.adj() += C_adj * AsolveB * 2 * arena_D.val_op().asDiagonal();
          arena_D.adj() += C_adj * BTAsolveB.diagonal();
        });

    return res;
  } else if (!is_constant<Ta>::value && !is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_D.asDiagonal() * arena_B.val_op().transpose() * AsolveB)
                  .trace();

    reverse_pass_callback([arena_A, AsolveB, arena_B, arena_D, res]() mutable {
      double C_adj = res.adj();

      arena_A.adj()
          -= C_adj * AsolveB * arena_D.asDiagonal() * AsolveB.transpose();
      arena_B.adj() += C_adj * AsolveB * 2 * arena_D.asDiagonal();
    });

    return res;
  } else if (!is_constant<Ta>::value && is_constant<Tb>::value
             && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(value_of(B_ref)));
    auto BTAsolveB = to_arena(value_of(B_ref).transpose() * AsolveB);

    var res = (arena_D.val().asDiagonal() * BTAsolveB).trace();

    reverse_pass_callback(
        [arena_A, BTAsolveB, AsolveB, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().asDiagonal()
                           * AsolveB.transpose();
          arena_D.adj() += C_adj * BTAsolveB.diagonal();
        });

    return res;
  } else if (!is_constant<Ta>::value && is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Ta>> arena_A = A.matrix();
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(value_of(B_ref)));

    var res = (arena_D.asDiagonal() * value_of(B_ref).transpose() * AsolveB)
                  .trace();

    reverse_pass_callback([arena_A, AsolveB, arena_D, res]() mutable {
      double C_adj = res.adj();

      arena_A.adj() -= C_adj * AsolveB * arena_D.val_op().asDiagonal()
                       * AsolveB.transpose();
    });

    return res;
  } else if (is_constant<Ta>::value && !is_constant<Tb>::value
             && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));
    auto BTAsolveB = to_arena(arena_B.val_op().transpose() * AsolveB);

    var res = (arena_D.val().asDiagonal() * BTAsolveB).trace();

    reverse_pass_callback(
        [BTAsolveB, AsolveB, arena_B, arena_D, res]() mutable {
          double C_adj = res.adj();

          arena_B.adj() += C_adj * AsolveB * 2 * arena_D.val_op().asDiagonal();
          arena_D.adj() += C_adj * BTAsolveB.diagonal();
        });

    return res;
  } else if (is_constant<Ta>::value && !is_constant<Tb>::value
             && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;
    arena_t<promote_scalar_t<double, Td>> arena_D = value_of(D);
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_D.asDiagonal() * arena_B.val_op().transpose() * AsolveB)
                  .trace();

    reverse_pass_callback([AsolveB, arena_B, arena_D, res]() mutable {
      arena_B.adj() += res.adj() * AsolveB * 2 * arena_D.asDiagonal();
    });

    return res;
  } else if (is_constant<Ta>::value && is_constant<Tb>::value
             && !is_constant<Td>::value) {
    const auto& B_ref = to_ref(B);
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    auto BTAsolveB = to_arena(value_of(B_ref).transpose()
                              * A.ldlt().solve(value_of(B_ref)));

    var res = (arena_D.val().asDiagonal() * BTAsolveB).trace();

    reverse_pass_callback([BTAsolveB, arena_D, res]() mutable {
      arena_D.adj() += res.adj() * BTAsolveB.diagonal();
    });

    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
