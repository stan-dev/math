#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsNegBinomial, derivatives_diff_sizes) {
  using stan::math::neg_binomial_lpmf;
  using stan::math::var;

  int N = 100;
  double mu_dbl = 1.5;
  std::vector<double> phi_dbl{2, 4, 6, 8};

  var mu(mu_dbl);
  std::vector<var> phi;
  for (double i : phi_dbl) {
    phi.push_back(var(i));
  }
  var val = neg_binomial_lpmf(N, mu, phi);

  std::vector<var> x{mu};
  std::vector<double> gradients;
  val.grad(x, gradients);

  double eps = 1e-6;
  double grad_diff = (neg_binomial_lpmf(N, mu_dbl + eps, phi_dbl)
                      - neg_binomial_lpmf(N, mu_dbl - eps, phi_dbl))
                     / (2 * eps);
  EXPECT_FLOAT_EQ(grad_diff, gradients[0]);
}
