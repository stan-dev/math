#ifndef STAN_MATH_PRIM_FUN_BETA_BINOMIAL_CDF_2_HPP
#define STAN_MATH_PRIM_FUN_BETA_BINOMIAL_CDF_2_HPP

#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <iostream>

namespace stan {
namespace math {

double beta_binomial_cdf_2(int n, int N, double a, double b) {
  using std::log;

  double lbeta_ab = lbeta(a, b);

  int Np1 = N + 1;
  double log_a = log(a);
  double log_bpn = log(b + N);
  double lgamma_ai = lgamma(a);
  double lgamma_bNi = lgamma(b + N);
  double lgamma_abN = lgamma(a + b + N);
  double log_Np1mi = 0.0;
  double log_i;

  Eigen::VectorXd log_cdf_i(n+1);

  for (int i = 0; i <= n; ++i) {
    if (i != 0) {
      log_i = log(i);
      log_Np1mi += log(Np1 - i) - log_i;
      lgamma_ai += log(a + i - 1);
      lgamma_bNi -= log(b + N - i);
    }

    log_cdf_i[i] = log_Np1mi + (((lgamma_ai + lgamma_bNi) - lgamma_abN) - lbeta_ab);
    std::cout << "i: " << i << ", " << exp(log_cdf_i[i]) << std::endl;
  }

  return exp(log_sum_exp(log_cdf_i));
}

}
}
#endif