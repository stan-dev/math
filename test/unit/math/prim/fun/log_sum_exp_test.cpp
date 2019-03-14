
#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
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

TEST(MathFunctions, log_sum_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_PRED1(boost::math::isnan<double>, stan::math::log_sum_exp(1.0, nan));

  EXPECT_PRED1(boost::math::isnan<double>, stan::math::log_sum_exp(nan, 1.0));

  EXPECT_PRED1(boost::math::isnan<double>, stan::math::log_sum_exp(nan, nan));
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

TEST(MathFunctions_mat, log_sum_exp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_sum_exp;

  Matrix<double, Dynamic, Dynamic> m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  test_log_sum_exp(m);

  Matrix<double, Dynamic, 1> v(3);
  v << 1, 2, 3;
  test_log_sum_exp(v);

  Matrix<double, Dynamic, 1> rv(3);
  rv << 1, 2, 3;
  test_log_sum_exp(rv);

  Matrix<double, Dynamic, Dynamic> m_trivial(1, 1);
  m_trivial << 2;
  EXPECT_FLOAT_EQ(2, log_sum_exp(m_trivial));
}
