#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
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
return_type_t<T_shape, T_inv_scale> neg_binomial_cdf(const T_n& n,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static constexpr const char* function = "neg_binomial_cdf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  if (size_zero(n, alpha, beta)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_alpha = stan::math::size(alpha);
  size_t size_n_alpha = max_size(n, alpha);
  size_t max_size_seq_view = max_size(n, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_alpha_vec(size_alpha);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_n,
                T_shape>
      digamma_sum_vec(size_n_alpha);

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; i++) {
      digamma_alpha_vec[i] = digamma(alpha_vec.val(i));
    }
    for (size_t i = 0; i < size_n_alpha; i++) {
      const T_partials_return n_dbl = n_vec.val(i);
      const T_partials_return alpha_dbl = alpha_vec.val(i);
      digamma_sum_vec[i] = digamma(n_dbl + alpha_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(1.0);
    }

    const T_partials_return n_dbl = n_vec.val(i);
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);
    const T_partials_return inv_beta_p1 = inv(beta_dbl + 1);
    const T_partials_return p_dbl = beta_dbl * inv_beta_p1;
    const T_partials_return d_dbl = square(inv_beta_p1);

    const T_partials_return P_i = inc_beta(alpha_dbl, n_dbl + 1.0, p_dbl);

    P *= P_i;

    if (!is_constant_all<T_shape>::value) {
      partials<0>(ops_partials)[i]
          += inc_beta_dda(alpha_dbl, n_dbl + 1, p_dbl, digamma_alpha_vec[i],
                          digamma_sum_vec[i])
             / P_i;
    }

    if (!is_constant_all<T_inv_scale>::value) {
      partials<1>(ops_partials)[i]
          += inc_beta_ddz(alpha_dbl, n_dbl + 1.0, p_dbl) * d_dbl / P_i;
    }
  }

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; ++i) {
      partials<0>(ops_partials)[i] *= P;
    }
  }

  if (!is_constant_all<T_inv_scale>::value) {
    for (size_t i = 0; i < stan::math::size(beta); ++i) {
      partials<1>(ops_partials)[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
