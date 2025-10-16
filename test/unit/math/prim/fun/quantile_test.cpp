#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

template <typename T>
void test_quantile_double() {
  using stan::math::index_type_t;
  using stan::math::quantile;

  // check size 0 argument throws
  T c0(0);
  EXPECT_NO_THROW(quantile(c0, 0.3));

  // check size 1 argument works
  T c(1);
  c[0] = 1.7;
  EXPECT_EQ(quantile(c, 0), 1.7);
  EXPECT_EQ(quantile(c, 1), 1.7);
  EXPECT_EQ(quantile(c, 0.33), 1.7);
  EXPECT_EQ(quantile(c, 0.68), 1.7);

  // check outside-of-[0,1] throws
  EXPECT_THROW(quantile(c, 1.1), std::domain_error);
  EXPECT_THROW(quantile(c, -0.1), std::domain_error);

  // check nan argument throws
  EXPECT_THROW(quantile(c, std::numeric_limits<double>::quiet_NaN()),
               std::domain_error);
  std::decay_t<T> qt_vec(5);
  stan::math::fill(qt_vec,
                   std::numeric_limits<stan::scalar_type_t<T>>::quiet_NaN());
  EXPECT_THROW(quantile(qt_vec, 0.5), std::domain_error);

  // check quantile results against R
  T v(5);
  v[0] = -0.07;
  v[1] = 0.35;
  v[2] = 2.65;
  v[3] = -0.28;
  v[4] = 1.23;
  EXPECT_FLOAT_EQ(quantile(v, 0), -0.28);
  EXPECT_FLOAT_EQ(quantile(v, 0.1), -0.196);
  EXPECT_FLOAT_EQ(quantile(v, 0.2), -0.112);
  EXPECT_FLOAT_EQ(quantile(v, 0.3), 0.014);
  EXPECT_FLOAT_EQ(quantile(v, 0.4), 0.182);
  EXPECT_FLOAT_EQ(quantile(v, 0.5), 0.350);
  EXPECT_FLOAT_EQ(quantile(v, 0.6), 0.702);
  EXPECT_FLOAT_EQ(quantile(v, 0.7), 1.054);
  EXPECT_FLOAT_EQ(quantile(v, 0.8), 1.514);
  EXPECT_FLOAT_EQ(quantile(v, 0.9), 2.082);
  EXPECT_FLOAT_EQ(quantile(v, 1), 2.65);
}

TEST(MathFunctions, quantileStdVecDouble) {
  test_quantile_double<std::vector<double>>();
}

TEST(MathFunctions, quantileEigenVectorXd) {
  test_quantile_double<Eigen::VectorXd>();
}

TEST(MathFunctions, quantileEigenRowVectorXd) {
  test_quantile_double<Eigen::RowVectorXd>();
}

TEST(MathFunctions, quantileStdVecInt) {
  using stan::math::quantile;

  std::vector<int> y0;
  EXPECT_NO_THROW(quantile(y0, 0.7));

  std::vector<int> v{-1, 3, 5, -10};

  EXPECT_FLOAT_EQ(quantile(v, 0), -10.0);
  EXPECT_FLOAT_EQ(quantile(v, 0.1), -7.3);
  EXPECT_FLOAT_EQ(quantile(v, 0.2), -4.6);
  EXPECT_FLOAT_EQ(quantile(v, 0.3), -1.9);
  EXPECT_FLOAT_EQ(quantile(v, 0.4), -0.2);
  EXPECT_FLOAT_EQ(quantile(v, 0.5), 1.0);
  EXPECT_FLOAT_EQ(quantile(v, 0.6), 2.2);
  EXPECT_FLOAT_EQ(quantile(v, 0.7), 3.2);
  EXPECT_FLOAT_EQ(quantile(v, 0.8), 3.8);
  EXPECT_FLOAT_EQ(quantile(v, 0.9), 4.4);
  EXPECT_FLOAT_EQ(quantile(v, 1), 5.0);
}

template <typename T, typename Tp, stan::require_all_vector_t<T, Tp>* = nullptr>
void test_quantile_double() {
  using stan::math::index_type_t;
  using stan::math::quantile;

  // check quantile results against R
  T v(5);
  v[0] = -0.07;
  v[1] = 0.35;
  v[2] = 2.65;
  v[3] = -0.28;
  v[4] = 1.23;

  Tp p(4);
  p[0] = 0;
  p[1] = 0.1;
  p[2] = 0.2;
  p[3] = 1;

  std::vector<double> ret = quantile(v, p);
  EXPECT_FLOAT_EQ(ret[0], -0.28);
  EXPECT_FLOAT_EQ(ret[1], -0.196);
  EXPECT_FLOAT_EQ(ret[2], -0.112);
  EXPECT_FLOAT_EQ(ret[3], 2.65);

  // check outside-of-[0,1] throws
  Tp p2(2);
  p2[0] = 0.5;
  p2[1] = 1.1;

  Tp p3(2);
  p3[0] = 0.5;
  p3[1] = -0.1;

  EXPECT_THROW(quantile(v, p2), std::domain_error);
  EXPECT_THROW(quantile(v, p3), std::domain_error);

  std::decay_t<T> qt_vec(5);
  stan::math::fill(qt_vec,
                   std::numeric_limits<stan::scalar_type_t<T>>::quiet_NaN());
  std::decay_t<Tp> pt_vec(4);
  stan::math::fill(pt_vec,
                   std::numeric_limits<stan::scalar_type_t<Tp>>::quiet_NaN());

  // check nan argument throws
  EXPECT_THROW(
      quantile(v, std::numeric_limits<stan::scalar_type_t<Tp>>::quiet_NaN()),
      std::domain_error);
  EXPECT_THROW(quantile(qt_vec, p), std::domain_error);

  // check size 0 in both arguments throws
  T v0(0);
  Tp p0(0);
  EXPECT_NO_THROW(quantile(v0, p));
  EXPECT_NO_THROW(quantile(v, p0));

  // check size 1 first argument works
  T v1(1);
  v1[0] = -0.07;
  std::vector<double> ret1 = quantile(v1, p);
  EXPECT_FLOAT_EQ(ret1[0], -0.07);
  EXPECT_FLOAT_EQ(ret1[1], -0.07);
  EXPECT_FLOAT_EQ(ret1[2], -0.07);
  EXPECT_FLOAT_EQ(ret1[3], -0.07);
}

TEST(MathFunctions, quantileStdVecDoubleStdVecDouble) {
  test_quantile_double<std::vector<double>, std::vector<double>>();
}

TEST(MathFunctions, quantileEigenVectorXdStdVecDouble) {
  test_quantile_double<Eigen::VectorXd, std::vector<double>>();
}

TEST(MathFunctions, quantileEigenRowVectorXdStdVecDouble) {
  test_quantile_double<Eigen::RowVectorXd, std::vector<double>>();
}
