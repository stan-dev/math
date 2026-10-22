#include <stan/math/fwd.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <cmath>

TEST(AgradFwd, log_softmax_large_dynamic_range) {
  using stan::math::fvar;
  using stan::math::log_softmax;
  using stan::math::vector_fd;

  Eigen::VectorXd val(3);
  val << 0, -100, -1000;

  Eigen::VectorXd d(3);
  d << 1, 2, 3;

  vector_fd x = stan::math::to_fvar(val, d);

  auto y = log_softmax(x);

  EXPECT_MATRIX_EQ(y.val(), val);
  EXPECT_MATRIX_EQ(y.d(), d.array() - (1.0 + 2.0 * std::exp(-100.0)));
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

TEST(AgradFwd, log_softmax_positive_infinity_consistent_with_prim) {
  using stan::math::fvar;
  using stan::math::log_softmax;

  constexpr double inf = std::numeric_limits<double>::infinity();

  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x(4);
  x(0) = fvar<double>(1.0, 1.0);
  x(1) = fvar<double>(2.0, 2.0);
  x(2) = fvar<double>(3.0, 3.0);
  x(3) = fvar<double>(inf, 4.0);

  auto y = log_softmax(x);

  EXPECT_EQ(-inf, y(0).val());
  EXPECT_EQ(-inf, y(1).val());
  EXPECT_EQ(-inf, y(2).val());
  EXPECT_TRUE(std::isnan(y(3).val()));
}
