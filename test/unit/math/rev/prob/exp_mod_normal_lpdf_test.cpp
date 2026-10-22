#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_propto_vectorized.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace {
auto exp_mod_normal_propto
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf<true>(y, mu, sigma, lambda);
      };

std::vector<stan::math::var> to_vars(const std::vector<double>& x) {
  return {x.begin(), x.end()};
}

// z = (y - mu) / sigma - lambda sigma spans the deep lower tail, the small
// branch, the positive tail and z > 40.
const std::vector<double> y{-40.0, -2.0, 0.5, 3.0, 10.0, 60.0};
const std::vector<double> mu{0.1, 0.3, -0.4, 1.0, 0.0, -0.2};
const std::vector<double> sigma{1.3, 0.5, 1.2, 2.0, 0.9, 1.0};
const std::vector<double> lambda{0.7, 4.0, 1.0, 0.2, 0.5, 2.5};
}  // namespace

TEST_F(AgradRev, ProbDistributionsExpModNormal_propto_vector_var) {
  using stan::test::expect_propto_vectorized;
  expect_propto_vectorized(exp_mod_normal_propto, to_vars(y), 0.1, 1.3, 0.7);
  expect_propto_vectorized(exp_mod_normal_propto, 0.5, to_vars(mu), 1.3, 0.7);
  expect_propto_vectorized(exp_mod_normal_propto, 0.5, 0.1, to_vars(sigma),
                           0.7);
  expect_propto_vectorized(exp_mod_normal_propto, 0.5, 0.1, 1.3,
                           to_vars(lambda));
  expect_propto_vectorized(exp_mod_normal_propto, to_vars(y), mu,
                           to_vars(sigma), 0.7);
  expect_propto_vectorized(exp_mod_normal_propto, y, to_vars(mu), sigma,
                           to_vars(lambda));
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_propto_scalar_var_broadcast) {
  using stan::math::var;
  using stan::test::expect_propto_vectorized;
  // A scalar var is shared by every element, so log(lambda) counts N times.
  expect_propto_vectorized(exp_mod_normal_propto, y, 0.1, 1.3, var(0.7));
  expect_propto_vectorized(exp_mod_normal_propto, y, mu, sigma, var(0.7));
  expect_propto_vectorized(exp_mod_normal_propto, var(0.5), 0.1, 1.3, lambda);
  expect_propto_vectorized(exp_mod_normal_propto, 0.5, var(0.1), 1.3, lambda);
  expect_propto_vectorized(exp_mod_normal_propto, y, mu, var(1.3), lambda);
}

TEST_F(AgradRev, ProbDistributionsExpModNormal_propto_eigen) {
  using stan::math::var;
  using stan::test::expect_propto_vectorized;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_d(y.size());
  Eigen::Matrix<var, 1, Eigen::Dynamic> lambda_v(lambda.size());
  for (size_t i = 0; i < y.size(); ++i) {
    y_d(i) = y[i];
    lambda_v(i) = lambda[i];
  }
  expect_propto_vectorized(exp_mod_normal_propto, y_d, var(0.1), 1.3, lambda_v);
}
