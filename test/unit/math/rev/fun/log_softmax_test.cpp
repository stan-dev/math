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

TEST(AgradRev, log_softmax_negative_infinity_has_finite_adjoint) {
  using stan::math::log_softmax;
  using stan::math::sum;
  using stan::math::var_value;
  
  const double neg_inf = -std::numeric_limits<double>::infinity();
  
  Eigen::VectorXd v(4);
  v << neg_inf, 1.0, 2.0, 3.0;
  
  var_value<Eigen::VectorXd> x(v);
  
  auto y = log_softmax(x);
  sum(y).grad();
  
  const double e1 = std::exp(1.0);
  const double e2 = std::exp(2.0);
  const double e3 = std::exp(3.0);
  const double denom = e1 + e2 + e3;
  
  const double p0 = 0.0;
  const double p1 = e1 / denom;
  const double p2 = e2 / denom;
  const double p3 = e3 / denom;
  
  EXPECT_NEAR(1.0, x.adj()(0), 1e-12);
  EXPECT_NEAR(1.0 - 4.0 * p1, x.adj()(1), 1e-12);
  EXPECT_NEAR(1.0 - 4.0 * p2, x.adj()(2), 1e-12);
  EXPECT_NEAR(1.0 - 4.0 * p3, x.adj()(3), 1e-12);
}

TEST(AgradRev, log_softmax_adjoint_uses_stable_softmax) {
  using stan::math::log_softmax;
  using stan::math::sum;
  using stan::math::var_value;
  
  Eigen::VectorXd v(4);
  
  v << 629.7901581243797,
       31.52411463,
       608.19720553,
       120.94829574;
  
  var_value<Eigen::VectorXd> x(v);
  auto y = log_softmax(x);
  
  sum(y).grad();
  
  // Independent max-shifted softmax calculation
  Eigen::VectorXd p = (v.array() -  v.maxCoeff()).exp();
  p /= p.sum();
  
  // d/dx_m sum_k log_softmax(x)_k
  // = 1 - n * softmax(x)_m.
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    const double expected = 1.0 - v.size() * p(i);
    
    EXPECT_NEAR(expected, x.adj()(i), 1e-12);
  }
}
