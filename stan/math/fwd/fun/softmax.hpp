#ifndef STAN_MATH_FWD_FUN_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_SOFTMAX_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/softmax.hpp>

namespace stan {
namespace math {

template <typename Container,
          require_vector_st<is_fvar, Container>...>
inline auto softmax(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [&](const auto& alpha) {
    using T_fvar_inner = typename value_type_t<decltype(alpha)>::Scalar;
    plain_type_t<decltype(alpha)> softmax_alpha(alpha.size());

    softmax_alpha.val() = softmax(alpha.val());
    softmax_alpha.d().setZero();

    for (int m = 0; m < alpha.size(); ++m) {
      T_fvar_inner negative_alpha_m_d_times_softmax_alpha_t_m
          = -alpha(m).d_ * softmax_alpha.val().coeff(m);
      for (int k = 0; k < alpha.size(); ++k) {
        if (m == k) {
          softmax_alpha(k).d_
              += softmax_alpha.val().coeff(k)
                 * (alpha(m).d_ + negative_alpha_m_d_times_softmax_alpha_t_m);
        } else {
          softmax_alpha(k).d_
              += negative_alpha_m_d_times_softmax_alpha_t_m
                  * softmax_alpha.val().coeff(k);
        }
      }
    }

    return softmax_alpha;
  });
}


}  // namespace math
}  // namespace stan
#endif
