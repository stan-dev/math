#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, squared_distance) {
  double x1 = 1;
  double x2 = 4;

  EXPECT_FLOAT_EQ(9, stan::math::squared_distance(x1, x2));
  EXPECT_FLOAT_EQ(9, stan::math::squared_distance(x2, x1));
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(x1, x1));
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(x2, x2));
}

TEST(MathFunctions, squared_distance_nan) {
  double x = 1;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::squared_distance(x, nan), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, x), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, nan), std::domain_error);
}

TEST(MathFunctions, squared_distance_inf) {
  double x = 1;
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::squared_distance(x, inf), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, x), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, inf), std::domain_error);
}

TEST(MathMatrixPrimMat, squared_distance_vector_vector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> v1, v2;

  v1.resize(3);
  v2.resize(3);
  v1 << 1, 3, -5;
  v2 << 4, -2, -1;

  EXPECT_FLOAT_EQ(50, stan::math::squared_distance(v1, v2));

  v1.resize(0);
  v2.resize(0);
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(v1, v2));

  v1.resize(1);
  v2.resize(2);
  v1 << 1;
  v2 << 2, 3;
  EXPECT_THROW(stan::math::squared_distance(v1, v2), std::invalid_argument);
}

TEST(MathMatrixPrimMat, squared_distance_rowvector_vector) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v;

  rv.resize(3);
  v.resize(3);
  rv << 1, 3, -5;
  v << 4, -2, -1;
  EXPECT_FLOAT_EQ(50, stan::math::squared_distance(rv, v));

  rv.resize(0);
  v.resize(0);
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(rv, v));

  rv.resize(1);
  v.resize(2);
  rv << 1;
  v << 2, 3;
  EXPECT_THROW(stan::math::squared_distance(rv, v), std::invalid_argument);
}

TEST(MathMatrixPrimMat, squared_distance_vector_rowvector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> v;
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv;

  v.resize(3);
  rv.resize(3);
  v << 1, 3, -5;
  rv << 4, -2, -1;
  EXPECT_FLOAT_EQ(50, stan::math::squared_distance(v, rv));

  v.resize(0);
  rv.resize(0);
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(v, rv));

  v.resize(1);
  rv.resize(2);
  v << 1;
  rv << 2, 3;
  EXPECT_THROW(stan::math::squared_distance(v, rv), std::invalid_argument);
}

TEST(MathMatrixPrimMat, squared_distance_special_values) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> v1, v2;
  v1.resize(1);
  v2.resize(1);

  v1 << 0;
  v2 << std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(stan::math::is_nan(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(stan::math::is_nan(stan::math::squared_distance(v2, v1)));

  v1 << 0;
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_inf(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(stan::math::is_inf(stan::math::squared_distance(v2, v1)));

  v1 << std::numeric_limits<double>::infinity();
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_nan(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(stan::math::is_nan(stan::math::squared_distance(v2, v1)));

  v1 << -std::numeric_limits<double>::infinity();
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(stan::math::is_inf(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(stan::math::is_inf(stan::math::squared_distance(v2, v1)));
}
