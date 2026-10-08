#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

template <typename T_shape>
inline return_type_t<double, T_shape> do_lkj_constant(const T_shape& eta,
                                                      const unsigned int& K) {
  // Lewandowski, Kurowicka, and Joe (2009) theorem 5
  using T_partials_return = partials_return_t<T_shape>;
  const T_partials_return eta_val = value_of(eta);
  T_partials_return constant;
  const int Km1 = K - 1;
  using stan::math::lgamma;
  if (eta_val == 1.0) {
    // C++ integer division is appropriate in this block
    Eigen::VectorXd denominator(Km1 / 2);
    for (int k = 1; k <= denominator.rows(); k++) {
      denominator(k - 1) = lgamma(2.0 * k);
    }
    constant = -denominator.sum();
    if ((K % 2) == 1) {
      constant -= 0.25 * (K * K - 1) * LOG_PI - 0.25 * (Km1 * Km1) * LOG_TWO
                  - Km1 * lgamma(0.5 * (K + 1));
    } else {
      constant -= 0.25 * K * (K - 2) * LOG_PI
                  + 0.25 * (3 * K * K - 4 * K) * LOG_TWO + K * lgamma(0.5 * K)
                  - Km1 * lgamma(static_cast<double>(K));
    }
  } else {
    // Km1 lgamma(eta + Km1 / 2) - sum_k lgamma(eta + (Km1 - k) / 2) is a sum
    // of differences lgamma(x + k / 2) - lgamma(x), x = eta + (Km1 - k) / 2,
    // which cancel for large eta. Each one equals
    // lgamma(k / 2) - lbeta(k / 2, x), which does not cancel.
    constant = -0.25 * Km1 * K * LOG_PI;
    for (int k = 1; k <= Km1; k++) {
      constant += lgamma(0.5 * k) - lbeta(0.5 * k, eta_val + 0.5 * (Km1 - k));
    }
  }
  // d/deta = sum_k digamma(x + k / 2) - digamma(x), also at eta = 1
  auto ops_partials = make_partials_propagator(eta);
  if constexpr (is_autodiff_v<T_shape>) {
    T_partials_return d_eta(0.0);
    for (int k = 1; k <= Km1; k++) {
      d_eta += digamma_diff(eta_val + 0.5 * (Km1 - k), 0.5 * k);
    }
    partials<0>(ops_partials)[0] = d_eta;
  }
  return ops_partials.build(constant);
}

// LKJ_Corr(y|eta) [ y correlation matrix (not covariance matrix)
//                  eta > 0; eta == 1 <-> uniform]
template <bool propto, typename T_y, typename T_shape>
inline return_type_t<T_y, T_shape> lkj_corr_lpdf(const T_y& y,
                                                 const T_shape& eta) {
  static constexpr const char* function = "lkj_corr_lpdf";

  return_type_t<T_y, T_shape> lp(0.0);
  const auto& y_ref = to_ref(y);
  check_positive(function, "Shape parameter", eta);
  check_corr_matrix(function, "Correlation matrix", y_ref);

  const unsigned int K = y.rows();
  if (K == 0) {
    return 0.0;
  }

  if constexpr (include_summand<propto, T_shape>::value) {
    lp += do_lkj_constant(eta, K);
  }
  if constexpr (is_constant_all<scalar_type<T_shape>>::value) {
    if (eta == 1.0) {
      return lp;
    }
  }

  if constexpr (!include_summand<propto, T_y, T_shape>::value) {
    return lp;
  }

  value_type_t<T_y> sum_values = sum(log(y_ref.ldlt().vectorD()));
  lp += (eta - 1.0) * sum_values;
  return lp;
}

template <typename T_y, typename T_shape>
inline return_type_t<T_y, T_shape> lkj_corr_lpdf(const T_y& y,
                                                 const T_shape& eta) {
  return lkj_corr_lpdf<false>(y, eta);
}

}  // namespace math
}  // namespace stan
#endif
