#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

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
          require_all_matrix_t<Ta, Tb>* = nullptr,
          require_return_type_t<is_var, Ta, Tb, Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<Ta> arena_a = a;
  arena_t<Tb> arena_b = b;
  return make_callback_var(
      hypergeometric_pFq(value_of(arena_a), value_of(arena_b), value_of(z)),
      [arena_a, arena_b, z](auto& vi) mutable {
        auto grad_tuple = grad_pFq(arena_a, arena_b, z);
        if (!is_constant<Ta>::value) {
          forward_as<promote_scalar_t<var, Ta>>(arena_a).adj()
              += vi.adj() * std::get<0>(grad_tuple);
        }
        if (!is_constant<Tb>::value) {
          forward_as<promote_scalar_t<var, Tb>>(arena_b).adj()
              += vi.adj() * std::get<1>(grad_tuple);
        }
        if (!is_constant<Tz>::value) {
          forward_as<promote_scalar_t<var, Tz>>(z).adj()
              += vi.adj() * std::get<2>(grad_tuple);
        }
      });
}
}  // namespace math
}  // namespace stan
#endif
