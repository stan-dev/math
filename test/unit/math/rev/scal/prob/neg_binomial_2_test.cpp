#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsNegBinomial, derivatives) {
  using stan::math::is_nan;
  using stan::math::neg_binomial_2_log;
  using stan::math::var;

  int N = 100;
  double mu_dbl = 8;
  double phi_dbl = 1.5;

  for (int k = 0; k < 20; ++k) {
    var mu(mu_dbl);
    var phi(phi_dbl);
    var val = neg_binomial_2_log(N, mu, phi);

    std::vector<var> x;
    x.push_back(mu);
    x.push_back(phi);

    std::vector<double> gradients;
    val.grad(x, gradients);

    for (int i = 0; i < 2; ++i) {
      EXPECT_FALSE(is_nan(gradients[i]));
    }

    std::vector<double> finite_diffs;
    double eps = 1e-10;
    double inv2e = 0.5 / eps;
    double dmu = neg_binomial_2_log(N, mu_dbl + eps, phi_dbl)
                 - neg_binomial_2_log(N, mu_dbl - eps, phi_dbl);
    double dphi = neg_binomial_2_log(N, mu_dbl, phi_dbl + eps)
                  - neg_binomial_2_log(N, mu_dbl, phi_dbl - eps);
    finite_diffs.push_back(dmu * inv2e);
    finite_diffs.push_back(dphi * inv2e);

    for (int i = 0; i < 2; ++i) {
      EXPECT_NEAR(gradients[i], finite_diffs[i], 1.0) << 
        "for mu = " << mu_dbl << "  +/- epsilon, phi = " << phi_dbl <<
        " +/- epsilon";
    }

    phi_dbl *= 10;
  }
}

TEST(ProbDistributionsNegativeBinomial2, proptoAtPoissonCutoff) {
  using stan::math::neg_binomial_2_lpmf;
  using stan::math::var;
  using stan::math::internal::neg_binomial_2_phi_cutoff;

  var mu_var(10);
  int y = 11.8;
  var value_before_cutoff = neg_binomial_2_lpmf<true, int, var, double>(
    y, mu_var, neg_binomial_2_phi_cutoff - 1e-8);
  var value_after_cutoff = neg_binomial_2_lpmf<true, int, var, double>(
    y, mu_var, neg_binomial_2_phi_cutoff + 1e-8);

  EXPECT_NEAR(value_of(value_before_cutoff), value_of(value_after_cutoff), 1);  
}
