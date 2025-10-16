#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsNegBinomial, derivatives_lcdf) {
  using stan::math::neg_binomial_lcdf;
  using stan::math::var;

  std::vector<double> N{1, 2, 3};
  double alpha_dbl = 8;
  double beta_dbl = 1.5;

  var alpha(alpha_dbl);
  var beta(beta_dbl);

  var val = neg_binomial_lcdf(N, alpha, beta);
  std::vector<var> x{alpha, beta};
  std::vector<double> gradients;
  val.grad(x, gradients);

  double epsilon = 1e-6;

  double grad_diff1 = (neg_binomial_lcdf(N, alpha_dbl + epsilon, beta_dbl)
                       - neg_binomial_lcdf(N, alpha_dbl - epsilon, beta_dbl))
                      / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff1, gradients[0]);

  double grad_diff2 = (neg_binomial_lcdf(N, alpha_dbl, beta_dbl + epsilon)
                       - neg_binomial_lcdf(N, alpha_dbl, beta_dbl - epsilon))
                      / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff2, gradients[1]);
}
