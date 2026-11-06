#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_propto_vectorized.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace {
auto skew_normal_propto
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf<true>(y, mu, sigma, alpha);
      };

std::vector<stan::math::var> to_vars(const std::vector<double>& x) {
  return {x.begin(), x.end()};
}

// alpha * z spans the deep lower tail, the small branch and z > 40.
const std::vector<double> y{-50.0, -2.0, 0.3, 1.5, 4.0, 50.0};
const std::vector<double> mu{0.1, -0.4, 0.6, 1.0, -1.0, 0.0};
const std::vector<double> sigma{1.0, 0.5, 1.2, 2.0, 0.9, 1.0};
const std::vector<double> alpha{1.0, -3.0, -0.5, 0.0, 2.5, -1.0};
}  // namespace

TEST_F(AgradRev, ProbDistributionsSkewNormal_propto_vector_var) {
  using stan::test::expect_propto_vectorized;
  expect_propto_vectorized(skew_normal_propto, to_vars(y), 0.2, 1.3, 0.7);
  expect_propto_vectorized(skew_normal_propto, 0.2, to_vars(mu), 1.3, 0.7);
  expect_propto_vectorized(skew_normal_propto, 0.2, 0.1, to_vars(sigma), 0.7);
  // z is a scalar, so only the broadcast log Phi(alpha z) terms remain.
  expect_propto_vectorized(skew_normal_propto, 0.2, 0.1, 1.3, to_vars(alpha));
  expect_propto_vectorized(skew_normal_propto, to_vars(y), mu, to_vars(sigma),
                           0.7);
  expect_propto_vectorized(skew_normal_propto, y, to_vars(mu), sigma,
                           to_vars(alpha));
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_propto_scalar_var_broadcast) {
  using stan::math::var;
  using stan::test::expect_propto_vectorized;
  // A scalar var is shared by every element, so its terms count N times.
  expect_propto_vectorized(skew_normal_propto, var(0.2), 0.1, 1.3, alpha);
  expect_propto_vectorized(skew_normal_propto, 0.2, var(0.1), 1.3, alpha);
  expect_propto_vectorized(skew_normal_propto, y, 0.1, var(1.3), 0.7);
  expect_propto_vectorized(skew_normal_propto, y, mu, var(1.3), alpha);
  expect_propto_vectorized(skew_normal_propto, y, mu, sigma, var(0.7));
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_propto_eigen) {
  using stan::math::var;
  using stan::test::expect_propto_vectorized;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v(y.size());
  Eigen::Matrix<double, 1, Eigen::Dynamic> alpha_d(alpha.size());
  for (size_t i = 0; i < y.size(); ++i) {
    y_v(i) = y[i];
    alpha_d(i) = alpha[i];
  }
  expect_propto_vectorized(skew_normal_propto, y_v, 0.1, var(1.3), alpha_d);
}
