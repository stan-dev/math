#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
// log Phi(x) through erfc, accurate for x >= -37.
double log_Phi_ref(double x) {
  return std::log(0.5 * std::erfc(-x / std::sqrt(2.0)));
}

double exp_mod_normal_ref(double y, double mu, double sigma, double lambda) {
  return std::log(lambda) + lambda * (mu - y + 0.5 * lambda * sigma * sigma)
         + log_Phi_ref((y - mu) / sigma - lambda * sigma);
}

void expect_near_rel(double expected, double actual) {
  EXPECT_NEAR(expected, actual, 1e-12 * std::max(1.0, std::fabs(expected)));
}

constexpr double mu = 0.1;
constexpr double sigma = 1.3;
constexpr double lambda = 0.7;

// Values of y whose z = (y - mu) / sigma - lambda sigma crosses
// +/-4 sqrt(2) and +/-0.46875 sqrt(2) and reaches the positive tail.
std::vector<double> branch_ys() {
  std::vector<double> y;
  for (double z : {-30.0, -5.7, -5.6, -0.67, 0.0, 0.66, 5.6, 5.7, 20.0, 30.0}) {
    y.push_back(mu + sigma * (z + lambda * sigma));
  }
  return y;
}
}  // namespace

TEST(ProbDistributionsExpModNormal, lpdf_matches_closed_form_across_branches) {
  using stan::math::exp_mod_normal_lpdf;
  const auto y = branch_ys();
  double sum_ref = 0;
  for (double y_i : y) {
    const double ref = exp_mod_normal_ref(y_i, mu, sigma, lambda);
    expect_near_rel(ref, exp_mod_normal_lpdf(y_i, mu, sigma, lambda));
    sum_ref += ref;
  }
  expect_near_rel(sum_ref, exp_mod_normal_lpdf(y, mu, sigma, lambda));
}

TEST(ProbDistributionsExpModNormal, lpdf_vector_lambda_broadcasts_scalars) {
  using stan::math::exp_mod_normal_lpdf;
  const std::vector<double> lambdas{0.05, 0.7, 2.5, 10.0, 20.0};
  double sum_ref = 0;
  for (double lambda_i : lambdas) {
    sum_ref += exp_mod_normal_ref(0.3, mu, sigma, lambda_i);
  }
  expect_near_rel(sum_ref, exp_mod_normal_lpdf(0.3, mu, sigma, lambdas));
}

TEST(ProbDistributionsExpModNormal, lpdf_vector_mu_sigma) {
  using stan::math::exp_mod_normal_lpdf;
  const std::vector<double> mus{-20.0, -1.0, 0.0, 2.0, 15.0};
  const std::vector<double> sigmas{0.5, 1.0, 2.0, 3.0, 0.8};
  double sum_mu = 0;
  double sum_both = 0;
  for (size_t i = 0; i < mus.size(); ++i) {
    sum_mu += exp_mod_normal_ref(0.3, mus[i], sigma, lambda);
    sum_both += exp_mod_normal_ref(0.3, mus[i], sigmas[i], lambda);
  }
  expect_near_rel(sum_mu, exp_mod_normal_lpdf(0.3, mus, sigma, lambda));
  expect_near_rel(sum_both, exp_mod_normal_lpdf(0.3, mus, sigmas, lambda));
}

TEST(ProbDistributionsExpModNormal, lpdf_size_one_vector_matches_scalar) {
  using stan::math::exp_mod_normal_lpdf;
  const std::vector<double> y{0.3};
  const Eigen::VectorXd lambda_vec = Eigen::VectorXd::Constant(1, lambda);
  EXPECT_FLOAT_EQ(exp_mod_normal_lpdf(0.3, mu, sigma, lambda),
                  exp_mod_normal_lpdf(y, mu, sigma, lambda));
  EXPECT_FLOAT_EQ(exp_mod_normal_lpdf(0.3, mu, sigma, lambda),
                  exp_mod_normal_lpdf(y, mu, sigma, lambda_vec));
}

TEST(ProbDistributionsExpModNormal, lpdf_propto_and_empty_return_zero) {
  using stan::math::exp_mod_normal_lpdf;
  const std::vector<double> y{-1.0, 0.3, 2.0};
  const std::vector<double> empty;
  EXPECT_EQ(0.0, exp_mod_normal_lpdf<true>(y, mu, sigma, lambda));
  EXPECT_EQ(0.0, exp_mod_normal_lpdf(empty, mu, sigma, lambda));
  EXPECT_EQ(0.0, exp_mod_normal_lpdf(empty, mu, sigma, empty));
}

TEST(ProbDistributionsExpModNormal, lpdf_vector_errors) {
  using stan::math::exp_mod_normal_lpdf;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> good{-1.0, 0.3, 2.0};

  EXPECT_THROW(exp_mod_normal_lpdf(std::vector<double>{-1.0, nan, 2.0}, mu,
                                   sigma, lambda),
               std::domain_error);
  for (double bad : {nan, inf, -inf}) {
    const std::vector<double> bad_mu{-1.0, bad, 2.0};
    EXPECT_THROW(exp_mod_normal_lpdf(good, bad_mu, sigma, lambda),
                 std::domain_error);
  }
  for (double bad : {nan, 0.0, -1.0, inf, -inf}) {
    const std::vector<double> bad_vec{0.5, bad, 2.0};
    EXPECT_THROW(exp_mod_normal_lpdf(good, mu, bad_vec, lambda),
                 std::domain_error);
    EXPECT_THROW(exp_mod_normal_lpdf(good, mu, sigma, bad_vec),
                 std::domain_error);
  }
  EXPECT_THROW(exp_mod_normal_lpdf(0.3, mu, inf, lambda), std::domain_error);
  // The checks run before the empty-input return.
  const std::vector<double> empty;
  EXPECT_THROW(exp_mod_normal_lpdf(empty, mu, sigma, -1.0), std::domain_error);
}

TEST(ProbDistributionsExpModNormal, lpdf_inconsistent_sizes) {
  using stan::math::exp_mod_normal_lpdf;
  const std::vector<double> three{-1.0, 0.3, 2.0};
  const std::vector<double> two{0.5, 1.0};
  EXPECT_THROW(exp_mod_normal_lpdf(two, three, sigma, lambda),
               std::invalid_argument);
  EXPECT_THROW(exp_mod_normal_lpdf(three, two, sigma, lambda),
               std::invalid_argument);
  EXPECT_THROW(exp_mod_normal_lpdf(three, mu, two, lambda),
               std::invalid_argument);
  EXPECT_THROW(exp_mod_normal_lpdf(three, mu, sigma, two),
               std::invalid_argument);
}
