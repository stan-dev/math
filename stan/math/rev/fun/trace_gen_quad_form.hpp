#ifndef STAN_MATH_REV_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/trace_gen_quad_form.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename Td, typename Ta, typename Tb,
          require_any_vt_var<Td, Ta, Tb>* = nullptr,
          require_all_eigen_t<Td, Ta, Tb>* = nullptr>
inline var trace_gen_quad_form(const Td& D, const Ta& A, const Tb& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

  using A_ref_t = ref_type_t<Ta>;
  using B_ref_t = ref_type_t<Tb>;
  using D_ref_t = ref_type_t<Td>;

  A_ref_t A_ref = A;
  B_ref_t B_ref = B;
  D_ref_t D_ref = D;

  arena_matrix<promote_scalar_t<double, Ta>> arena_A_val;
  arena_matrix<promote_scalar_t<double, Tb>> arena_B_val = value_of(B_ref);
  arena_matrix<promote_scalar_t<double, Td>> arena_D_val;

  arena_matrix<promote_scalar_t<var, Ta>> arena_A;
  arena_matrix<promote_scalar_t<var, Tb>> arena_B;
  arena_matrix<promote_scalar_t<var, Td>> arena_D;

  if (!is_constant<Ta>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<Tb>::value) {
    arena_B = B_ref;
  }

  if (!is_constant<Td>::value) {
    arena_D = D_ref;
  }

  if (!is_constant_all<Tb, Td>::value) {
    arena_A_val = value_of(A_ref);
  }

  if (!is_constant_all<Ta, Tb>::value) {
    arena_D_val = value_of(D_ref);
  }

  var res;

  if (is_constant<Ta>::value && is_constant<Tb>::value
      && !is_constant<Td>::value) {
    res = (value_of(D) * arena_B_val.transpose() * arena_A_val * arena_B_val)
              .trace();
  } else if (!is_constant<Ta>::value && is_constant<Tb>::value
             && is_constant<Td>::value) {
    res = (arena_D_val * arena_B_val.transpose() * value_of(A) * arena_B_val)
              .trace();
  } else {
    res = (arena_D_val * arena_B_val.transpose() * arena_A_val * arena_B_val)
              .trace();
  }

  reverse_pass_callback([arena_A, arena_B, arena_D, arena_A_val, arena_B_val,
                         arena_D_val, res]() mutable {
    double C_adj = res.adj();

    if (!is_constant<Ta>::value)
      arena_A.adj() += C_adj * arena_B_val * arena_D_val.transpose()
                       * arena_B_val.transpose();

    if (!is_constant<Tb>::value)
      arena_B.adj() += C_adj
                       * (arena_A_val * arena_B_val * arena_D_val
                          + arena_A_val.transpose() * arena_B_val
                                * arena_D_val.transpose());

    if (!is_constant<Td>::value)
      arena_D.adj()
          += C_adj * ((arena_A_val * arena_B_val).transpose() * arena_B_val);
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
