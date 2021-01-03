#ifndef STAN_MATH_PRIM_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/trace.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <exception>

namespace stan {
namespace math {

/**
 * Return the trace of D times the quadratic form of B and A.
 * That is, `trace_gen_quad_form(D, A, B) = trace(D * B' * A * B).`
 *
 * @tparam TD type of the first matrix or expression
 * @tparam TA type of the second matrix or expression
 * @tparam TB type of the third matrix or expression
 *
 * @param D multiplier
 * @param A outside term in quadratic form
 * @param B inner term in quadratic form
 * @return trace(D * B' * A * B)
 * @throw std::domain_error if A or D is not square
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by D.
 */
template <typename TD, typename TA, typename TB,
          typename = require_all_eigen_t<TD, TA, TB>,
          typename = require_all_not_vt_var<TD, TA, TB>,
          typename = require_any_not_vt_arithmetic<TD, TA, TB>>
inline auto trace_gen_quad_form(const TD& D, const TA& A, const TB& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);
  const auto& B_ref = to_ref(B);
  return multiply(B_ref, D.transpose()).cwiseProduct(multiply(A, B_ref)).sum();
}

/**
 * Return the trace of D times the quadratic form of B and A.
 * That is, `trace_gen_quad_form(D, A, B) = trace(D * B' * A * B).`
 * This is the overload for arithmetic types to allow Eigen's expression
 * templates to be used for efficiency.
 *
 * @tparam EigMatD type of the first matrix or expression
 * @tparam EigMatA type of the second matrix or expression
 * @tparam EigMatB type of the third matrix or expression
 *
 * @param D multiplier
 * @param A outside term in quadratic form
 * @param B inner term in quadratic form
 * @return trace(D * B' * A * B)
 * @throw std::domain_error if A or D is not square
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by D.
 */
template <typename EigMatD, typename EigMatA, typename EigMatB,
          require_all_eigen_vt<std::is_arithmetic, EigMatD, EigMatA,
                               EigMatB>* = nullptr>
inline double trace_gen_quad_form(const EigMatD& D, const EigMatA& A,
                                  const EigMatB& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);
  const auto& B_ref = to_ref(B);
  return (B_ref * D.transpose()).cwiseProduct(A * B_ref).sum();
}

}  // namespace math
}  // namespace stan

#endif
