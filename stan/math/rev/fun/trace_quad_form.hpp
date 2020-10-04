#ifndef STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/trace_quad_form.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline var trace_quad_form(const EigMat1& A, const EigMat2& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  auto arena_B_val = to_arena(value_of(B_ref));
  var res;
  if (!is_constant<EigMat1>::value && !is_constant<EigMat2>::value) {
    arena_t<promote_scalar_t<var, EigMat1>> arena_A = A_ref;
    arena_t<promote_scalar_t<var, EigMat2>> arena_B = B_ref;
    auto arena_A_val = to_arena(value_of(arena_A));

    res = (arena_B_val.transpose() * arena_A_val * arena_B_val).trace();

    reverse_pass_callback([arena_A, arena_B, arena_A_val, arena_B_val,
                           res]() mutable {
      arena_A.adj() += res.adj() * arena_B_val * arena_B_val.transpose();
      arena_B.adj()
          += res.adj() * (arena_A_val + arena_A_val.transpose()) * arena_B_val;
    });
  } else if (!is_constant<EigMat1>::value) {
    arena_t<promote_scalar_t<var, EigMat1>> arena_A = A_ref;

    res = (arena_B_val.transpose() * value_of(A_ref) * arena_B_val).trace();

    reverse_pass_callback([arena_A, arena_B_val, res]() mutable {
      arena_A.adj() += res.adj() * arena_B_val * arena_B_val.transpose();
    });
  } else {
    arena_t<promote_scalar_t<var, EigMat2>> arena_B = B_ref;
    auto arena_A_val = to_arena(value_of(A_ref));

    res = (arena_B_val.transpose() * arena_A_val * arena_B_val).trace();

    reverse_pass_callback([arena_B, arena_A_val, arena_B_val, res]() mutable {
      double C_adj = res.adj();

      arena_B.adj()
          += res.adj() * (arena_A_val + arena_A_val.transpose()) * arena_B_val;
    });
  }

  return res;
}

}  // namespace math
}  // namespace stan
#endif
