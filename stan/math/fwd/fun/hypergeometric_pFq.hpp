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
 * @tparam FvarVec1 An Eigen vector with arithmetic scalar type
 * @tparam FvarVec2 An Eigen vector with arithmetic scalar type
 * @tparam Fvar An arithmetic scalar
 * @param[in] a Vector of 'a' arguments (of length p)
 * @param[in] b Vector of 'b' arguments (of length q)
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
template <typename FvarVec1, typename FvarVec2, typename Fvar,
          require_all_matrix_t<FvarVec1, FvarVec2>* = nullptr,
          require_return_type_t<is_fvar, FvarVec1, FvarVec2, Fvar>* = nullptr>
inline return_type_t<FvarVec1, FvarVec2, Fvar> hypergeometric_pFq(const FvarVec1& a, const FvarVec2& b,
                                                    const Fvar& z) {
  using fvar_t = return_type_t<FvarVec1, FvarVec2, Fvar>;
  ref_type_t<FvarVec1> a_ref = a;
  ref_type_t<FvarVec2> b_ref = b;
  auto grad_tuple = grad_pFq(a_ref, b_ref, z);

  typename fvar_t::Scalar grad = 0;

  if (!is_constant<FvarVec1>::value) {
    grad += dot_product(forward_as<promote_scalar_t<fvar_t, FvarVec1>>(a_ref).d(),
                        std::get<0>(grad_tuple));
  }
  if (!is_constant<FvarVec2>::value) {
    grad += dot_product(forward_as<promote_scalar_t<fvar_t, FvarVec2>>(b_ref).d(),
                        std::get<1>(grad_tuple));
  }
  if (!is_constant<Fvar>::value) {
    grad += forward_as<promote_scalar_t<fvar_t, Fvar>>(z).d_
            * std::get<2>(grad_tuple);
  }

  return fvar_t(
      hypergeometric_pFq(value_of(a_ref), value_of(b_ref), value_of(z)), grad);
}

}  // namespace math
}  // namespace stan
#endif
