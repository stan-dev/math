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
 * @tparam EigMat1 type of the first matrix
 * @tparam T2 type of elements in the LDLT_factor
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam EigMat3 type of the second matrix
 *
 * @param D a square matrix
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <typename Td, typename Ta, typename Tb, int R, int C,
          require_all_eigen_t<Td, Tb>* = nullptr,
          require_any_st_var<Td, Ta, Tb>* = nullptr>
inline var trace_gen_inv_quad_form_ldlt(const Td& D,
                                        const LDLT_factor<Ta, R, C>& A,
                                        const Tb& B) {
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A, "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  if (D.size() == 0 || A.cols() == 0) {
    return 0;
  }

  const auto& B_ref = to_ref(B);
  const auto& D_ref = to_ref(D);

  auto arena_B_val = to_arena(value_of(B_ref));
  arena_t<promote_scalar_t<double, Td>> arena_D_val;
  auto AsolveB = to_arena(A.solve(arena_B_val));
  auto BTAsolveB = to_arena(arena_B_val.transpose() * AsolveB);

  arena_t<promote_scalar_t<var, Tb>> arena_B;
  arena_t<promote_scalar_t<var, Td>> arena_D;

  if (!is_constant<Tb>::value) {
    arena_B = B_ref;
  }

  if (!is_constant<Td>::value) {
    arena_D = D_ref;
  }

  if (!is_constant_all<Ta, Tb>::value) {
    arena_D_val = value_of(D_ref);
  }

  var res;

  if (!is_constant_all<Ta, Tb>::value) {
    res = (arena_D_val * BTAsolveB).trace();
  } else {
    res = (value_of(D) * BTAsolveB).trace();
  }

  reverse_pass_callback([A, BTAsolveB, AsolveB, arena_B, arena_D, arena_B_val, arena_D_val,
                         res]() mutable {
    double C_adj = res.adj();

    if (!is_constant<Ta>::value) {
      forward_as<LDLT_factor<var, R, C>>(A).alloc_->arena_A_.adj()
          += -C_adj * AsolveB * arena_D_val.transpose() * AsolveB.transpose();
    }

    if (!is_constant<Tb>::value)
      arena_B.adj()
          += C_adj * AsolveB * (arena_D_val + arena_D_val.transpose());

    if (!is_constant<Td>::value)
      arena_D.adj() += C_adj * BTAsolveB;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
