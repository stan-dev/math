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

double skew_normal_ref(double y, double mu, double sigma, double alpha) {
  const double z = (y - mu) / sigma;
  return stan::math::LOG_TWO - stan::math::HALF_LOG_TWO_PI - std::log(sigma)
         - 0.5 * z * z + log_Phi_ref(alpha * z);
}

void expect_near_rel(double expected, double actual) {
  EXPECT_NEAR(expected, actual, 1e-12 * std::max(1.0, std::fabs(expected)));
}

// alpha * z for alpha = +/-1 crosses +/-4 sqrt(2) and +/-0.46875 sqrt(2).
const std::vector<double> branch_points{-30.0, -5.7, -5.6, -0.67, 0.0,
                                        0.66,  5.6,  5.7,  20.0,  30.0};
}  // namespace

TEST(ProbDistributionsSkewNormal, lpdf_matches_closed_form_across_branches) {
  using stan::math::skew_normal_lpdf;
  for (double alpha : {-1.0, 1.0}) {
    double sum_ref = 0;
    for (double y : branch_points) {
      const double ref = skew_normal_ref(y, 0.0, 1.0, alpha);
      expect_near_rel(ref, skew_normal_lpdf(y, 0.0, 1.0, alpha));
      sum_ref += ref;
    }
    expect_near_rel(sum_ref, skew_normal_lpdf(branch_points, 0.0, 1.0, alpha));
  }
}

TEST(ProbDistributionsSkewNormal, lpdf_vector_alpha_broadcasts_scalars) {
  using stan::math::skew_normal_lpdf;
  // z = 1 is shared, so its quadratic and log(sigma) terms count N times.
  double sum_ref = 0;
  for (double alpha : branch_points) {
    sum_ref += skew_normal_ref(1.3, 0.3, 1.0, alpha);
  }
  expect_near_rel(sum_ref, skew_normal_lpdf(1.3, 0.3, 1.0, branch_points));
}

TEST(ProbDistributionsSkewNormal, lpdf_zero_shape_is_normal) {
  using stan::math::normal_lpdf;
  using stan::math::skew_normal_lpdf;
  const std::vector<double> y{-2.0, 0.3, 4.0};
  const std::vector<double> sigma{1e-8, 1.0, 1e8};
  expect_near_rel(normal_lpdf(y, 0.1, sigma),
                  skew_normal_lpdf(y, 0.1, sigma, 0.0));
}

TEST(ProbDistributionsSkewNormal, lpdf_large_shape_is_twice_normal) {
  using stan::math::LOG_TWO;
  using stan::math::normal_lpdf;
  using stan::math::skew_normal_lpdf;
  // alpha * z > 40, where Phi(alpha z) is 1.
  const std::vector<double> y{1.0, 2.0, 3.0};
  expect_near_rel(normal_lpdf(y, 0.0, 1.0) + 3 * LOG_TWO,
                  skew_normal_lpdf(y, 0.0, 1.0, 100.0));
}

TEST(ProbDistributionsSkewNormal, lpdf_size_one_vector_matches_scalar) {
  using stan::math::skew_normal_lpdf;
  const std::vector<double> y{0.3};
  const Eigen::VectorXd alpha = Eigen::VectorXd::Constant(1, 0.5);
  EXPECT_FLOAT_EQ(skew_normal_lpdf(0.3, 0.1, 1.2, 0.5),
                  skew_normal_lpdf(y, 0.1, 1.2, 0.5));
  EXPECT_FLOAT_EQ(skew_normal_lpdf(0.3, 0.1, 1.2, 0.5),
                  skew_normal_lpdf(y, 0.1, 1.2, alpha));
}

TEST(ProbDistributionsSkewNormal, lpdf_propto_and_empty_return_zero) {
  using stan::math::skew_normal_lpdf;
  const std::vector<double> y{-1.0, 0.3, 2.0};
  const std::vector<double> empty;
  EXPECT_EQ(0.0, skew_normal_lpdf<true>(y, 0.1, 1.2, 0.5));
  EXPECT_EQ(0.0, skew_normal_lpdf(empty, 0.1, 1.2, 0.5));
  EXPECT_EQ(0.0, skew_normal_lpdf(empty, 0.1, 1.2, empty));
}

TEST(ProbDistributionsSkewNormal, lpdf_vector_errors) {
  using stan::math::skew_normal_lpdf;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> good{-1.0, 0.3, 2.0};
  const std::vector<double> good_sigma{0.5, 1.0, 2.0};

  EXPECT_THROW(
      skew_normal_lpdf(std::vector<double>{-1.0, nan, 2.0}, 0.1, 1.2, 0.5),
      std::domain_error);
  for (double bad : {nan, inf, -inf}) {
    const std::vector<double> bad_vec{-1.0, bad, 2.0};
    EXPECT_THROW(skew_normal_lpdf(good, bad_vec, 1.2, 0.5), std::domain_error);
    EXPECT_THROW(skew_normal_lpdf(good, 0.1, 1.2, bad_vec), std::domain_error);
  }
  for (double bad : {nan, 0.0, -1.0, -inf}) {
    const std::vector<double> bad_sigma{0.5, bad, 2.0};
    EXPECT_THROW(skew_normal_lpdf(good, 0.1, bad_sigma, 0.5),
                 std::domain_error);
  }
  // The checks run before the empty-input return.
  const std::vector<double> empty;
  EXPECT_THROW(skew_normal_lpdf(empty, 0.1, -1.0, 0.5), std::domain_error);
}

TEST(ProbDistributionsSkewNormal, lpdf_inconsistent_sizes) {
  using stan::math::skew_normal_lpdf;
  const std::vector<double> three{-1.0, 0.3, 2.0};
  const std::vector<double> two{0.5, 1.0};
  EXPECT_THROW(skew_normal_lpdf(two, three, 1.2, 0.5), std::invalid_argument);
  EXPECT_THROW(skew_normal_lpdf(three, two, 1.2, 0.5), std::invalid_argument);
  EXPECT_THROW(skew_normal_lpdf(three, 0.1, two, 0.5), std::invalid_argument);
  EXPECT_THROW(skew_normal_lpdf(three, 0.1, 1.2, two), std::invalid_argument);
}
