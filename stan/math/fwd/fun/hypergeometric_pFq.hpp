#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalized hypergeometric (pFq) function applied to the
 * input arguments.
 *
 * @tparam Ta Type of Eigen vector with scalar type fvar or arithmetic
 * @tparam Tb Type of Eigen vector with scalar type fvar or arithmetic
 * @tparam Tz Scalar of type fvar or arithmetic
 * @param[in] a Vector of 'a' arguments (of length p)
 * @param[in] b Vector of 'b' arguments (of length q)
 * @param[in] z Scalar z argument
 * @return Generalized hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          typename FvarT = return_type_t<Ta, Tb, Tz>,
          bool grad_a = !is_constant<Ta>::value,
          bool grad_b = !is_constant<Tb>::value,
          bool grad_z = !is_constant<Tz>::value,
          require_all_vector_t<Ta, Tb>* = nullptr,
          require_fvar_t<FvarT>* = nullptr>
inline FvarT hypergeometric_pFq(Ta&& a, Tb&& b, Tz&& z) {
  auto&& a_ref = to_ref(as_column_vector_or_scalar(a));
  auto&& b_ref = to_ref(as_column_vector_or_scalar(b));
  auto&& a_val = value_of(a_ref);
  auto&& b_val = value_of(b_ref);
  auto&& z_val = value_of(z);

  partials_type_t<FvarT> pfq_val = hypergeometric_pFq(a_val, b_val, z_val);
  auto grad_tuple
      = grad_pFq<grad_a, grad_b, grad_z>(pfq_val, a_val, b_val, z_val);

  FvarT rtn = FvarT(pfq_val, 0.0);

  if constexpr (grad_a) {
    rtn.d_ += dot_product(a_ref.d(), std::get<0>(grad_tuple));
  }
  if constexpr (grad_b) {
    rtn.d_ += dot_product(b_ref.d(), std::get<1>(grad_tuple));
  }
  if constexpr (grad_z) {
    rtn.d_ += z.d_ * std::get<2>(grad_tuple);
  }

  return rtn;
}

}  // namespace math
}  // namespace stan
#endif
