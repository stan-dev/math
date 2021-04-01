#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> neg_binomial_lcdf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta_param) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static const char* function = "neg_binomial_lcdf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta_param);
  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta_param;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  if (size_zero(n, alpha, beta_param)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_alpha_ref, T_beta_ref> ops_partials(alpha_ref,
                                                              beta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_n = stan::math::size(n);
  size_t size_alpha = stan::math::size(alpha);
  size_t size_n_alpha = max_size(n, alpha);
  size_t max_size_seq_view = max_size(n, alpha, beta_param);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_n; i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(negative_infinity());
    }
  }

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_n>
      digammaN_vec(size_n);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digammaAlpha_vec(size_alpha);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_n,
                T_shape>
      digammaSum_vec(size_n_alpha);

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_n; i++) {
      digammaN_vec[i] = digamma(n_vec.val(i) + 1);
    }
    for (size_t i = 0; i < size_alpha; i++) {
      digammaAlpha_vec[i] = digamma(alpha_vec.val(i));
    }
    for (size_t i = 0; i < size_n_alpha; i++) {
      const T_partials_return n_dbl = n_vec.val(i);
      const T_partials_return alpha_dbl = alpha_vec.val(i);
      digammaSum_vec[i] = digamma(n_dbl + alpha_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(0.0);
    }

    const T_partials_return n_dbl = n_vec.val(i);
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);
    const T_partials_return inv_beta_p1 = inv(beta_dbl + 1);
    const T_partials_return p_dbl = beta_dbl * inv_beta_p1;
    const T_partials_return d_dbl = square(inv_beta_p1);
    const T_partials_return Pi = inc_beta(alpha_dbl, n_dbl + 1.0, p_dbl);
    const T_partials_return beta_func = beta(n_dbl + 1, alpha_dbl);

    P += log(Pi);

    if (!is_constant_all<T_shape>::value) {
      T_partials_return g1 = 0;
      T_partials_return g2 = 0;

      grad_reg_inc_beta(g1, g2, alpha_dbl, n_dbl + 1, p_dbl,
                        digammaAlpha_vec[i], digammaN_vec[i], digammaSum_vec[i],
                        beta_func);
      ops_partials.edge1_.partials_[i] += g1 / Pi;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[i] += d_dbl * pow(1 - p_dbl, n_dbl)
                                          * pow(p_dbl, alpha_dbl - 1)
                                          / (beta_func * Pi);
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
