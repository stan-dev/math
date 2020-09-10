#ifndef STAN_MATH_REV_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Compute the trace of an inverse quadratic form premultiplied by a
 * square matrix. This computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 *
 * @tparam T type of elements in the LDLT_factor
 * @tparam R number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam EigMat type of the first matrix
 *
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <int R, int C, typename T, typename EigMat,
          require_eigen_t<EigMat>* = nullptr,
          require_any_st_var<T, EigMat>* = nullptr>
inline var trace_inv_quad_form_ldlt(const LDLT_factor<T, R, C>& A,
                                    const EigMat& B) {
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  if (A.rows() == 0)
    return 0.0;

  using B_ref_t = ref_type_t<EigMat>;

  B_ref_t B_ref = B;
  arena_matrix<promote_scalar_t<double, EigMat>> arena_B_val = value_of(B_ref);
  arena_matrix<promote_scalar_t<var, EigMat>> arena_B;
  arena_matrix<Eigen::Matrix<double, R, EigMat::ColsAtCompileTime>> AsolveB
      = A.solve(arena_B_val);

  if (!is_constant<EigMat>::value) {
    arena_B = B_ref;
  }

  var res = (arena_B_val.transpose() * AsolveB).trace();

  reverse_pass_callback([A, AsolveB, arena_B, arena_B_val, res]() mutable {
    double C_adj = res.adj();

    if (!is_constant<T>::value) {
      forward_as<const LDLT_factor<var, R, C>>(A).alloc_->arena_A_.adj()
          += -C_adj * AsolveB * AsolveB.transpose();
    }

    if (!is_constant<EigMat>::value)
      arena_B.adj() += 2 * C_adj * AsolveB;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
