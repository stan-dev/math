#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsNegBinomial2, derivatives) {
  using stan::math::is_nan;
  using stan::math::neg_binomial_2_lpmf;
  using stan::math::var;

  int N = 100;
  double mu_dbl = 8;
  double phi_dbl = 1.5;

  for (int k = 0; k < 20; ++k) {
    var mu(mu_dbl);
    var phi(phi_dbl);
    var val = neg_binomial_2_lpmf(N, mu, phi);

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

    double dmu = neg_binomial_2_lpmf(N, mu_dbl + eps, phi_dbl)
                 - neg_binomial_2_lpmf(N, mu_dbl - eps, phi_dbl);
    double dphi = neg_binomial_2_lpmf(N, mu_dbl, phi_dbl + eps)
                  - neg_binomial_2_lpmf(N, mu_dbl, phi_dbl - eps);
    finite_diffs.push_back(dmu * inv2e);
    finite_diffs.push_back(dphi * inv2e);

    for (int i = 0; i < 2; ++i) {
      EXPECT_NEAR(gradients[i], finite_diffs[i], 1.0);
    }

    phi_dbl *= 10;
  }
}

TEST(ProbDistributionsNegBinomial2, derivatives_diff_sizes) {
  using stan::math::neg_binomial_2_lpmf;
  using stan::math::var;

  int N = 100;
  double mu_dbl = 1.5;
  std::vector<double> phi_dbl{2, 4, 6, 8};

  var mu(mu_dbl);
  std::vector<var> phi;
  for (double i : phi_dbl) {
    phi.push_back(var(i));
  }
  var val = neg_binomial_2_lpmf(N, mu, phi);

  std::vector<var> x{mu};
  std::vector<double> gradients;
  val.grad(x, gradients);

  double eps = 1e-6;
  double grad_diff = (neg_binomial_2_lpmf(N, mu_dbl + eps, phi_dbl)
                      - neg_binomial_2_lpmf(N, mu_dbl - eps, phi_dbl))
                     / (2 * eps);
  EXPECT_FLOAT_EQ(grad_diff, gradients[0]);
}
