#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

#include <cmath>

TEST(AgradFwd, log_softmax_large_dynamic_range) {
  using stan::math::fvar;
  using stan::math::log_softmax;
  
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x(3);
  
  x(0) = fvar<double>(0.0, 1.0);
  x(1) = fvar<double>(-100.0, 2.0);
  x(2) = fvar<double>(-1000.0, 3.0);
  
  auto y = log_softmax(x);
  
  /*
   * log_softmax(x) = x - log_sum_exp(x).
   *
   * exp(-1000) underflows to zero in double precision. Therefore:
   *
   *   log(softmax(x)[2])
   *
   * gives -inf while
   *
   *   x[2] - log_sum_exp(x)
   *
   * correctly retains the log-probability.
   */
  
  EXPECT_TRUE(std::isfinite(y.val()(0)));
  EXPECT_TRUE(std::isfinite(y.val()(1)));
  EXPECT_TRUE(std::isfinite(y.val()(2)));
  
  EXPECT_DOUBLE_EQ(0.0, y.val()(0));
  EXPECT_DOUBLE_EQ(-100.0, y.val()(1));
  EXPECT_DOUBLE_EQ(-1000.0, y.val()(2));
  
  const double p1 = std::exp(-100.0);
  const double dot_sd = 1.0 + 2.0 * p1;
  
  EXPECT_NEAR(1.0 - dot_sd, y(0).d(), 1e-12);
  EXPECT_NEAR(2.0 - dot_sd, y(1).d(), 1e-12);
  EXPECT_NEAR(3.0 - dot_sd, y(2).d(), 1e-12);
}

TEST(AgradFwd, log_softmax_many_small_probabilities) {
  using stan::math::fvar;
  using stan::math::log_softmax;
  
  constexpr int n = 501;
  
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x(n);
  
  x(0) = fvar<double>(0.0, 1.0);
  
  for (int i = 1; i < n; ++i) {
    x(i) = fvar<double>(-36.0, 2.0);
  }
  
  const auto y = log_softmax(x);
  
  const double lse = std::log1p(500.0 * std::exp(-36.0));
  EXPECT_NEAR(-lse, y(0).val(), 1e-14);
  EXPECT_NEAR(-36.0 - lse, y(1).val(), 1e-14);
  EXPECT_NEAR(-36.0 - lse, y(500).val(), 1e-14);
  
  for (int i = 0; i < n; ++i) {
    EXPECT_TRUE(std::isfinite(y(i).val()));
    EXPECT_TRUE(std::isfinite(y(i).d()));
  }
  
  const double p_small = std::exp(-36.0) / (1.0 + 500.0 * std::exp(-36.0));
  const double p_large = 1.0 / (1.0 + 500.0 * std::exp(-36.0));
  
  const double tangent = p_large + 500.0 * 2.0 * p_small;
  
  EXPECT_NEAR(1.0 - tangent, y(0).d(), 1e-12);
  EXPECT_NEAR(2.0 - tangent, y(1).d(), 1e-12);
}
TEST(Agrad, log_softmax_positive_infinity_consistent_with_prim) {
  using stan::math::fvar;
  using stan::math::log_softmax;
  
  constexpr double inf = std::numeric_limits<double>::infinity();
  
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x(4);
  x(0) = fvar<double>(1.0, 1.0);
  x(1) = fvar<double>(2.0, 2.0);
  x(2) = fvar<double>(3.0, 3.0);
  x(3) = fvar<double>(inf, 4.0);
  
  auto y = log_softmax(x);
  
  // The value is defined as
  //
  //   x - log_sum_exp(x).
  //
  // Since log_sum_exp([1, 2, 3, +Inf]) = +Inf:
  //
  //   1 - Inf = -Inf
  //   2 - Inf = -Inf
  //   3 - Inf = -Inf
  //   Inf - Inf = NaN
  //
  
  EXPECT_EQ(-inf, y(0).val());
  EXPECT_EQ(-inf, y(1).val());
  EXPECT_EQ(-inf, y(2).val());
  EXPECT_TRUE(std::isnan(y(3).val()));
}