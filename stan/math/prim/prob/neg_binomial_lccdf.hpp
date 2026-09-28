#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_shape, typename T_inv_scale>
inline return_type_t<T_shape, T_inv_scale> neg_binomial_lccdf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta_param) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static constexpr const char* function = "neg_binomial_lccdf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta_param);

  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta_param;

  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  if (size_zero(n, alpha, beta_param)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_n = stan::math::size(n);
  size_t max_size_seq_view = max_size(n, alpha, beta_param);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_n; i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(negative_infinity());
    }

    const auto n_dbl = n_vec.val(i);
    const auto alpha_dbl = alpha_vec.val(i);
    const auto beta_dbl = beta_vec.val(i);
    const auto inv_beta_p1 = inv(beta_dbl + 1);
    const auto p_dbl = beta_dbl * inv_beta_p1;
    const auto d_dbl = square(inv_beta_p1);
    // The complement by the symmetry relation I_p(a, b) = I_{1-p}(b, a),
    // with 1 - p = 1 / (beta + 1) taken from inv_beta_p1 rather than as
    // 1 - p_dbl. Forming the complement as 1 - I_p(a, b) loses every digit
    // once I_p(a, b) rounds to 1, which happens for every complement below
    // eps.
    const auto Pi = inc_beta(n_dbl + 1.0, alpha_dbl, inv_beta_p1);
    const auto beta_func = beta(n_dbl + 1, alpha_dbl);

    P += log(Pi);

    // On the reflected arguments alpha is the second shape of inc_beta, and
    // only its partial is needed; it is the derivative of Pi itself.
    if constexpr (is_autodiff_v<T_shape>) {
      partials<0>(ops_partials)[i]
          += grad_reg_inc_beta<GRAD_RIB::G2>(
                 n_dbl + 1, alpha_dbl, inv_beta_p1, 0.0, digamma(alpha_dbl),
                 digamma(n_dbl + alpha_dbl + 1), 0.0)
             / Pi;
    }
    if constexpr (is_autodiff_v<T_inv_scale>) {
      // inc_beta_ddz is the same density, evaluated without the product of
      // two powers over a beta function, each of which can underflow on
      // its own.
      partials<1>(ops_partials)[i]
          -= d_dbl * inc_beta_ddz(alpha_dbl, n_dbl + 1.0, p_dbl) / Pi;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
