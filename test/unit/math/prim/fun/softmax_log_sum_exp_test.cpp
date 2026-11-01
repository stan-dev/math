#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, softmax_log_sum_exp_matches_separate_calls) {
  using stan::math::log_sum_exp;
  using stan::math::softmax;
  using stan::math::internal::softmax_log_sum_exp;

  Eigen::VectorXd v(4);
  v << 629.7901581243797, 31.52411463, 608.19720553, 120.94829574;
  auto res = softmax_log_sum_exp(v);
  EXPECT_MATRIX_EQ(softmax(v), res.first);
  EXPECT_FLOAT_EQ(log_sum_exp(v), res.second);

  Eigen::RowVectorXd rv(3);
  rv << -1.0, 1.0, 10.0;
  auto res_rv = softmax_log_sum_exp(rv);
  EXPECT_MATRIX_EQ(softmax(rv), res_rv.first);
  EXPECT_FLOAT_EQ(log_sum_exp(rv), res_rv.second);

  Eigen::MatrixXd m(2, 2);
  m << 0.0, -100.0, -1000.0, 3.0;
  auto res_m = softmax_log_sum_exp(m);
  EXPECT_MATRIX_EQ(softmax(m), res_m.first);
  EXPECT_FLOAT_EQ(log_sum_exp(m), res_m.second);
}

TEST(MathFunctions, softmax_log_sum_exp_empty) {
  using stan::math::internal::softmax_log_sum_exp;

  Eigen::VectorXd v(0);
  auto res = softmax_log_sum_exp(v);
  EXPECT_EQ(0, res.first.size());
  EXPECT_EQ(stan::math::NEGATIVE_INFTY, res.second);
}

TEST(MathFunctions, softmax_log_sum_exp_non_finite) {
  using stan::math::internal::softmax_log_sum_exp;
  constexpr double inf = std::numeric_limits<double>::infinity();

  Eigen::VectorXd v(3);
  v << 1.0, -inf, 2.0;
  auto res = softmax_log_sum_exp(v);
  EXPECT_FLOAT_EQ(2.0 + std::log(std::exp(-1.0) + 1.0), res.second);
  EXPECT_EQ(0.0, res.first(1));

  v << 1.0, inf, 2.0;
  res = softmax_log_sum_exp(v);
  EXPECT_EQ(inf, res.second);
  EXPECT_TRUE(res.first.array().isNaN().all());

  v << -inf, -inf, -inf;
  res = softmax_log_sum_exp(v);
  EXPECT_EQ(-inf, res.second);
  EXPECT_TRUE(res.first.array().isNaN().all());
}
