#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

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
  using stan::math::log_sum_exp;
  using Eigen::Dynamic;
  using Eigen::Matrix;

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
  //test_log_sum_exp(m);

  std::vector<Eigen::MatrixXd> st_m{m1,m2};
  std::vector<std::vector<double>> st_st{st1,st2};
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
