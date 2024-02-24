#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

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
template <
    typename Ta, typename Tb, typename Tz,
    typename FvarT = return_type_t<Ta, Tb, Tz>,
    bool GradA = !is_constant<Ta>::value,
    bool GradB = !is_constant<Tb>::value,
    bool GradZ = !is_constant<Tz>::value,
    require_all_vector_t<Ta, Tb>* = nullptr,
    require_return_type_t<is_fvar, FvarT>* = nullptr>
inline FvarT hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  using PartialsT = partials_type_t<FvarT>;
  using ARefT = ref_type_t<Ta>;
  using BRefT = ref_type_t<Tb>;

  ARefT a_ref = a;
  BRefT b_ref = b;
  auto&& a_val = value_of(a_ref);
  auto&& b_val = value_of(b_ref);
  auto&& z_val = value_of(z);
  PartialsT pfq_val = hypergeometric_pFq(a_val, b_val, z_val);
  auto grad_tuple = grad_pFq<GradA, GradB, GradZ>(pfq_val, a_val, b_val, z_val);

  FvarT rtn = FvarT(pfq_val, 0.0);

  if (GradA) {
    rtn.d_ += dot_product(forward_as<promote_scalar_t<FvarT, ARefT>>(a_ref).d(),
                            std::get<0>(grad_tuple));
  }
  if (GradB) {
    rtn.d_ += dot_product(forward_as<promote_scalar_t<FvarT, BRefT>>(b_ref).d(),
                            std::get<1>(grad_tuple));
  }
  if (GradZ) {
    rtn.d_ += forward_as<promote_scalar_t<FvarT, Tz>>(z).d_
                * std::get<2>(grad_tuple);
  }

  return rtn;
}

}  // namespace math
}  // namespace stan
#endif
