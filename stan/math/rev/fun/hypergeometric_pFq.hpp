#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalized hypergeometric function (pFq) applied to the
 * input arguments.
 *
 * @tparam Ta Type of Eigen vector with scalar type var or arithmetic
 * @tparam Tb Type of Eigen vector with scalar type var or arithmetic
 * @tparam Tz Scalar of type var or arithmetic
 * @param[in] a Vector of 'a' arguments (of length p)
 * @param[in] b Vector of 'b' arguments (of length q)
 * @param[in] z Scalar z argument
 * @return Generalized hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          bool grad_a = !is_constant<Ta>::value,
          bool grad_b = !is_constant<Tb>::value,
          bool grad_z = !is_constant<Tz>::value,
          require_all_vector_t<Ta, Tb>* = nullptr,
          require_return_type_t<is_var, Ta, Tb, Tz>* = nullptr>
inline var hypergeometric_pFq(Ta&& a, Tb&& b, Tz&& z) {
  auto&& arena_a = to_arena(as_column_vector_or_scalar(std::forward<Ta>(a)));
  auto&& arena_b = to_arena(as_column_vector_or_scalar(std::forward<Tb>(b)));
  auto pfq_val = hypergeometric_pFq(arena_a.val(), arena_b.val(), value_of(z));
  return make_callback_var(
      pfq_val, [arena_a, arena_b, z, pfq_val](auto& vi) mutable {
        auto grad_tuple = grad_pFq<grad_a, grad_b, grad_z>(
            pfq_val, arena_a.val(), arena_b.val(), value_of(z));
        if constexpr (grad_a) {
          arena_a.adj() += vi.adj() * std::get<0>(grad_tuple);
        }
        if constexpr (grad_b) {
          arena_b.adj() += vi.adj() * std::get<1>(grad_tuple);
        }
        if constexpr (grad_z) {
          z.adj() += vi.adj() * std::get<2>(grad_tuple);
        }
      });
}
}  // namespace math
}  // namespace stan
#endif
