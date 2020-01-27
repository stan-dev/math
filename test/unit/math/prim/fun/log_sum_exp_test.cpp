#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

void test_log_sum_exp(double a, double b) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;
  EXPECT_FLOAT_EQ(log(exp(a) + exp(b)), log_sum_exp(a, b));
}

void test_log_sum_exp(const std::vector<double>& as) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;
  double sum_exp = 0.0;
  for (size_t n = 0; n < as.size(); ++n)
    sum_exp += exp(as[n]);
  EXPECT_FLOAT_EQ(log(sum_exp), log_sum_exp(as));
}

TEST(MathFunctions, log_sum_exp) {
  using stan::math::log_sum_exp;
  std::vector<double> as;
  test_log_sum_exp(as);
  as.push_back(0.0);
  test_log_sum_exp(as);
  as.push_back(1.0);
  test_log_sum_exp(as);
  as.push_back(-1.0);
  test_log_sum_exp(as);
  as.push_back(-10000.0);
  test_log_sum_exp(as);

  as.push_back(10000.0);
  EXPECT_FLOAT_EQ(10000.0, log_sum_exp(as));
}

TEST(MathFunctions, log_sum_exp_2) {
  using stan::math::log_sum_exp;
  test_log_sum_exp(1.0, 2.0);
  test_log_sum_exp(1.0, 1.0);
  test_log_sum_exp(3.0, 2.0);
  test_log_sum_exp(-20.0, 12);
  test_log_sum_exp(-20.0, 12);

  // exp(10000.0) overflows
  EXPECT_FLOAT_EQ(10000.0, log_sum_exp(10000.0, 0.0));
  EXPECT_FLOAT_EQ(0.0, log_sum_exp(-10000.0, 0.0));
}

TEST(MathFunctions, log_sum_exp_2_inf) {
  using stan::math::log_sum_exp;
  double inf = std::numeric_limits<double>::infinity();
  test_log_sum_exp(1.0, -inf);
  test_log_sum_exp(-inf, 3.0);
  test_log_sum_exp(-inf, -inf);
  EXPECT_FLOAT_EQ(inf, log_sum_exp(inf, 3.0));
  EXPECT_FLOAT_EQ(inf, log_sum_exp(inf, inf));
}

TEST(MathFunctions, log_sum_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log_sum_exp(1.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::log_sum_exp(nan, 1.0)));

  EXPECT_TRUE(std::isnan(stan::math::log_sum_exp(nan, nan)));
}

template <int R, int C>
void test_log_sum_exp(const Eigen::Matrix<double, R, C>& as) {
  using stan::math::log_sum_exp;
  using std::exp;
  using std::log;
  double sum_exp = 0.0;
  for (int n = 0; n < as.size(); ++n)
    sum_exp += exp(as(n));
  EXPECT_FLOAT_EQ(log(sum_exp), log_sum_exp(as));
}

TEST(MathFunctions, log_sum_exp_mat) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_sum_exp;

  Matrix<double, Dynamic, Dynamic> m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;
  std::vector<double> st1{1, 2, 3, 4, 5, 6};
  Matrix<double, Dynamic, Dynamic> m2(3, 2);
  m2 << -1, -2, -3, -4, -5, -6;
  std::vector<double> st2{-1, -2, -3, -4, -5, -6};
  double m1_out = log_sum_exp(m1);
  double m2_out = log_sum_exp(m2);
  double st1_out = log_sum_exp(st1);
  double st2_out = log_sum_exp(st2);
  EXPECT_FLOAT_EQ(m1_out, st1_out);
  EXPECT_FLOAT_EQ(m2_out, st2_out);

  std::vector<Eigen::MatrixXd> st_m{m1, m2};
  std::vector<std::vector<double>> st_st{st1, st2};
  std::vector<double> m_out = log_sum_exp(st_m);
  std::vector<double> st_out = log_sum_exp(st_st);
  EXPECT_FLOAT_EQ(m1_out, m_out[0]);
  EXPECT_FLOAT_EQ(m2_out, m_out[1]);
  EXPECT_FLOAT_EQ(m1_out, st_out[0]);
  EXPECT_FLOAT_EQ(m2_out, st_out[1]);

  Matrix<double, Dynamic, 1> v(3);
  v << 1, 2, 3;
  test_log_sum_exp(v);

  Matrix<double, 1, Dynamic> rv(3);
  rv << 1, 2, 3;
  test_log_sum_exp(rv);

  Matrix<double, Dynamic, Dynamic> m_trivial(1, 1);
  m_trivial << 2;
  EXPECT_FLOAT_EQ(2, log_sum_exp(m_trivial));

  Matrix<double, Dynamic, 1> i(3);
  i << 1, 2, -std::numeric_limits<double>::infinity();
  test_log_sum_exp(i);

  Matrix<double, Dynamic, 1> ii(1);
  ii << -std::numeric_limits<double>::infinity();
  test_log_sum_exp(ii);
}
