#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

// Direct exercise of the var_value vector overload in
// stan/math/rev/fun/log_softmax.hpp for both row and column var_value vectors.

TEST_F(AgradRev, log_softmax_var_value_col_vector) {
  using stan::math::log_softmax;
  using stan::math::sum;
  using stan::math::var_value;

  Eigen::VectorXd v(3);
  v << -1, 1, 2;

  var_value<Eigen::VectorXd> x(v);
  auto y = log_softmax(x);

  EXPECT_EQ(3, y.val().rows());
  EXPECT_EQ(1, y.val().cols());

  const double lse = std::log(std::exp(-1.0) + std::exp(1.0) + std::exp(2.0));
  EXPECT_FLOAT_EQ(-1.0 - lse, y.val()(0));
  EXPECT_FLOAT_EQ(1.0 - lse, y.val()(1));
  EXPECT_FLOAT_EQ(2.0 - lse, y.val()(2));

  // d/dx[m] sum(y) = sum_k(delta_km - softmax[m]) = 1 - n * softmax[m]
  sum(y).grad();
  const double denom = std::exp(-1.0) + std::exp(1.0) + std::exp(2.0);
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(-1.0) / denom, x.adj()(0));
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(1.0) / denom, x.adj()(1));
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(2.0) / denom, x.adj()(2));
}

TEST_F(AgradRev, log_softmax_var_value_row_vector) {
  using stan::math::log_softmax;
  using stan::math::sum;
  using stan::math::var_value;

  Eigen::RowVectorXd v(3);
  v << -1, 1, 2;

  var_value<Eigen::RowVectorXd> x(v);
  auto y = log_softmax(x);

  // Output should preserve row shape.
  EXPECT_EQ(1, y.val().rows());
  EXPECT_EQ(3, y.val().cols());

  const double lse = std::log(std::exp(-1.0) + std::exp(1.0) + std::exp(2.0));
  EXPECT_FLOAT_EQ(-1.0 - lse, y.val()(0));
  EXPECT_FLOAT_EQ(1.0 - lse, y.val()(1));
  EXPECT_FLOAT_EQ(2.0 - lse, y.val()(2));

  // Same gradient formula; one softmax over the entire row vector.
  sum(y).grad();
  const double denom = std::exp(-1.0) + std::exp(1.0) + std::exp(2.0);
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(-1.0) / denom, x.adj()(0));
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(1.0) / denom, x.adj()(1));
  EXPECT_FLOAT_EQ(1.0 - 3.0 * std::exp(2.0) / denom, x.adj()(2));
}
