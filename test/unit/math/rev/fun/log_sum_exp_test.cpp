#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, log_sum_exp_tests_large_values) {
  using stan::math::var;

  // check autodiffing works with var types with large values
  var a = 1e50;
  var output = stan::math::log_sum_exp(a, a);
  output.grad();
  EXPECT_FLOAT_EQ(output.val(), log(2.0) + value_of(a));
  EXPECT_FLOAT_EQ(a.adj(), 1.0);

  var a2 = 1;
  var a3 = 1e50;
  var output2 = stan::math::log_sum_exp(a2, a3);
  output2.grad();
  EXPECT_FLOAT_EQ(a2.adj(), 0.0);
  EXPECT_FLOAT_EQ(a3.adj(), 1.0);

  var a4 = 1e50;
  var a5 = 1;
  var output3 = stan::math::log_sum_exp(a4, a5);
  output3.grad();
  EXPECT_FLOAT_EQ(a4.adj(), 1.0);
  EXPECT_FLOAT_EQ(a5.adj(), 0.0);

  // check autodiffing works with var types with large values
  var b = 1e20;
  var output6 = stan::math::log_sum_exp(b, b);
  output6.grad();
  EXPECT_FLOAT_EQ(output6.val(), log(2.0) + value_of(b));
  EXPECT_FLOAT_EQ(b.adj(), 1.0);

  var b2 = -2;
  var b3 = 1e20;
  var output7 = stan::math::log_sum_exp(b2, b3);
  output7.grad();
  EXPECT_FLOAT_EQ(b2.adj(), 0.0);
  EXPECT_FLOAT_EQ(b3.adj(), 1.0);

  var b4 = 1e20;
  var b5 = -2;
  var output8 = stan::math::log_sum_exp(b4, b5);
  output8.grad();
  EXPECT_FLOAT_EQ(b4.adj(), 1.0);
  EXPECT_FLOAT_EQ(b5.adj(), 0.0);

  // check arguement combinations of vars and doubles
  var a6 = 1e50;
  double a7 = 1;
  var output4 = stan::math::log_sum_exp(a6, a7);
  output4.grad();
  EXPECT_FLOAT_EQ(a6.adj(), 1.0);

  var a8 = 1;
  double a9 = 1e50;
  var output5 = stan::math::log_sum_exp(a8, a9);
  output5.grad();
  EXPECT_FLOAT_EQ(a8.adj(), 0.0);
}

TEST_F(AgradRev, log_sum_exp_negative_infinity_has_zero_adjoint) {
  using stan::math::log_sum_exp;
  using stan::math::var_value;
  
  const double neg_inf = -std::numeric_limits<double>::infinity();
  
  Eigen::VectorXd v(4);
  v << neg_inf, 1.0, 2.0, 3.0;
  
  var_value<Eigen::VectorXd> x(v);
  
  auto y = log_sum_exp(x);
  y.grad();
  
  const double e1 = std::exp(1.0);
  const double e2 = std::exp(2.0);
  const double e3 = std::exp(3.0);
  const double denom = e1 + e2 + e3;
  
  EXPECT_DOUBLE_EQ(0.0, x.adj()(0));
  EXPECT_NEAR(e1 / denom, x.adj()(1), 1e-12);
  EXPECT_NEAR(e2 / denom, x.adj()(2), 1e-12);
  EXPECT_NEAR(e3 / denom, x.adj()(3), 1e-12);
}

TEST_F(AgradRev, log_sum_exp_adjoint_uses_stable_softmax) {
  using stan::math::log_sum_exp;
  using stan::math::var_value;
  
  Eigen::VectorXd v(4);
  
  v << 629.7901581243797,
       31.52411463,
       608.19720553,
       120.94829574;
  
  var_value<Eigen::VectorXd> x(v);
  
  auto y = log_sum_exp(x);
  y.grad();
  
  // Independent max-shifted softmax calculation
  Eigen::VectorXd p = (v.array() - v.maxCoeff()).exp();
  p /= p.sum();
  
  // The gradient of log_sum_exp is softmax.
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    EXPECT_NEAR(p(i), x.adj()(i), 1e-12);
  }
}
