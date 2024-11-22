#ifndef STAN_MATH_REV_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/trace_inv_quad_form_ldlt.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Compute the trace of an inverse quadratic form premultiplied by a
 * square matrix. This computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 *
 * @tparam T1 type of elements in the LDLT_factor
 * @tparam T2 type of the second matrix
 *
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline var trace_inv_quad_form_ldlt(LDLT_factor<T1>& A, const T2& B) {
  check_multiplicable("trace_quad_form", "A", A.matrix(), "B", B);

  if (A.matrix().size() == 0)
    return 0.0;

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A.matrix();
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_B.val_op().transpose() * AsolveB).trace();

    reverse_pass_callback([arena_A, AsolveB, arena_B, res]() mutable {
      arena_A.adj() += -res.adj() * AsolveB * AsolveB.transpose();
      arena_B.adj() += 2 * res.adj() * AsolveB;
    });

    return res;
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A.matrix();
    const auto& B_ref = to_ref(B);

    auto AsolveB = to_arena(A.ldlt().solve(value_of(B_ref)));

    var res = (value_of(B_ref).transpose() * AsolveB).trace();

    reverse_pass_callback([arena_A, AsolveB, res]() mutable {
      arena_A.adj() += -res.adj() * AsolveB * AsolveB.transpose();
    });

    return res;
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    auto AsolveB = to_arena(A.ldlt().solve(arena_B.val()));

    var res = (arena_B.val_op().transpose() * AsolveB).trace();

    reverse_pass_callback([AsolveB, arena_B, res]() mutable {
      arena_B.adj() += 2 * res.adj() * AsolveB;
    });

    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
