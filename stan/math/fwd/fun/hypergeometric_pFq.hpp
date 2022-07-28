#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric (pFq) function applied to the
 * input arguments.
 *
 * @tparam Ta Type of Eigen vector with scalar type fvar or arithmetic
 * @tparam Tb Type of Eigen vector with scalar type fvar or arithmetic
 * @tparam Tz Scalar of type fvar or arithmetic
 * @param[in] a Vector of 'a' arguments (of length p)
 * @param[in] b Vector of 'b' arguments (of length q)
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_matrix_t<Ta, Tb>* = nullptr,
          require_return_type_t<is_fvar, Ta, Tb, Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  ref_type_t<Ta> a_ref = a;
  ref_type_t<Tb> b_ref = b;
  auto grad_tuple = grad_pFq(a_ref, b_ref, z);

  typename fvar_t::Scalar grad = 0;

  if (!is_constant<Ta>::value) {
    grad += dot_product(forward_as<promote_scalar_t<fvar_t, Ta>>(a_ref).d(),
                        std::get<0>(grad_tuple));
  }
  if (!is_constant<Tb>::value) {
    grad += dot_product(forward_as<promote_scalar_t<fvar_t, Tb>>(b_ref).d(),
                        std::get<1>(grad_tuple));
  }
  if (!is_constant<Tz>::value) {
    grad += forward_as<promote_scalar_t<fvar_t, Tz>>(z).d_
            * std::get<2>(grad_tuple);
  }

  return fvar_t(
      hypergeometric_pFq(value_of(a_ref), value_of(b_ref), value_of(z)), grad);
}

}  // namespace math
}  // namespace stan
#endif
