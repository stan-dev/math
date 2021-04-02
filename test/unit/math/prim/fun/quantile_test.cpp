#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

template <typename T>
void test_quantile_double() {
  using stan::math::index_type_t;
  using stan::math::quantile;

  T c(1);
  c[0] = 1.7;
  EXPECT_EQ(quantile(c, 0), 1.7);
  EXPECT_EQ(quantile(c, 1), 1.7);
  EXPECT_EQ(quantile(c, 0.33), 1.7);
  EXPECT_EQ(quantile(c, 0.68), 1.7);

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
  test_quantile_double<std::vector<double> >();
}

TEST(MathFunctions, quantileEigenVectorXd) {
  test_quantile_double<Eigen::VectorXd>();
}

TEST(MathFunctions, quantileEigenRowVectorXd) {
  test_quantile_double<Eigen::RowVectorXd>();
}