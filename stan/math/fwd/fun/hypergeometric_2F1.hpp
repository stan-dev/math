#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_2F1_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_2F1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/hypergeometric_2F1.hpp>
#include <stan/math/prim/fun/grad_2F1.hpp>

namespace stan {
namespace math {

/**
 * Returns the Gauss hypergeometric function applied to the
 * input arguments:
 * \f$_2F_1(a_1,a_2;b;z)\f$
 *
 * See 'grad_2F1.hpp' for the derivatives wrt each parameter
 *
 * @tparam Ta1 Type of scalar first 'a' argument
 * @tparam Ta2 Type of scalar second 'a' argument
 * @tparam Tb Type of scalar 'b' argument
 * @tparam Tz Type of scalar 'z' argument
 * @param[in] a1 First of 'a' arguments to function
 * @param[in] a2 Second of 'a' arguments to function
 * @param[in] b 'b' argument to function
 * @param[in] z Scalar z argument
 * @return Gauss hypergeometric function
 */
template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz>* = nullptr,
          require_any_fvar_t<Ta1, Ta2, Tb, Tz>* = nullptr>
inline return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1& a1,
                                                          const Ta2& a2,
                                                          const Tb& b,
                                                          const Tz& z) {
  using fvar_t = return_type_t<Ta1, Ta1, Tb, Tz>;

  auto a1_val = value_of(a1);
  auto a2_val = value_of(a2);
  auto b_val = value_of(b);
  auto z_val = value_of(z);

  auto grad_tuple = grad_2F1(a1, a2, b, z);

  typename fvar_t::Scalar grad = 0;

  if (!is_constant<Ta1>::value) {
    grad += forward_as<fvar_t>(a1).d() * std::get<0>(grad_tuple);
  }
  if (!is_constant<Ta2>::value) {
    grad += forward_as<fvar_t>(a2).d() * std::get<1>(grad_tuple);
  }
  if (!is_constant<Tb>::value) {
    grad += forward_as<fvar_t>(b).d() * std::get<2>(grad_tuple);
  }
  if (!is_constant<Tz>::value) {
    grad += forward_as<fvar_t>(z).d() * std::get<3>(grad_tuple);
  }

  return fvar_t(hypergeometric_2F1(a1_val, a2_val, b_val, z_val), grad);
}

}  // namespace math
}  // namespace stan
#endif
