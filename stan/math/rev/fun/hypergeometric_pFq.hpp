#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function (pFq) applied to the
 * input arguments.
 *
 * @tparam VarVec1 An Eigen vector with arithmetic scalar type
 * @tparam VarVec2 An Eigen vector with arithmetic scalar type
 * @tparam Var An arithmetic scalar
 * @param[in] a Vector of 'a' arguments (of length p)
 * @param[in] b Vector of 'b' arguments (of length q)
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
template <typename VarVec1, typename VarVec2, typename Var,
          require_all_matrix_t<VarVec1, VarVec2>* = nullptr,
          require_return_type_t<is_var, VarVec1, VarVec2, Var>* = nullptr>
inline var hypergeometric_pFq(const VarVec1& a, const VarVec2& b,
                              const Var& z) {
  arena_t<VarVec1> arena_a = a;
  arena_t<VarVec2> arena_b = b;
  return make_callback_var(
      hypergeometric_pFq(value_of(arena_a), value_of(arena_b), value_of(z)),
      [arena_a, arena_b, z](auto& vi) mutable {
        auto grad_tuple = grad_pFq(arena_a, arena_b, z);
        if (!is_constant<VarVec1>::value) {
          forward_as<promote_scalar_t<var, VarVec1>>(arena_a).adj()
              += vi.adj() * std::get<0>(grad_tuple);
        }
        if (!is_constant<VarVec2>::value) {
          forward_as<promote_scalar_t<var, VarVec2>>(arena_b).adj()
              += vi.adj() * std::get<1>(grad_tuple);
        }
        if (!is_constant<Var>::value) {
          forward_as<promote_scalar_t<var, Var>>(z).adj()
              += vi.adj() * std::get<2>(grad_tuple);
        }
      });
}
}  // namespace math
}  // namespace stan
#endif
