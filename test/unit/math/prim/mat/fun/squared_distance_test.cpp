#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <limits>

TEST(MathMatrix, squared_distance_vector_vector) {
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

TEST(MathMatrix, squared_distance_rowvector_vector) {
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

TEST(MathMatrix, squared_distance_vector_rowvector) {
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

TEST(MathMatrix, squared_distance_special_values) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> v1, v2;
  v1.resize(1);
  v2.resize(1);

  v1 << 0;
  v2 << std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(boost::math::isnan(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(boost::math::isnan(stan::math::squared_distance(v2, v1)));

  v1 << 0;
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(boost::math::isinf(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(boost::math::isinf(stan::math::squared_distance(v2, v1)));

  v1 << std::numeric_limits<double>::infinity();
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(boost::math::isnan(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(boost::math::isnan(stan::math::squared_distance(v2, v1)));

  v1 << -std::numeric_limits<double>::infinity();
  v2 << std::numeric_limits<double>::infinity();
  EXPECT_TRUE(boost::math::isinf(stan::math::squared_distance(v1, v2)));
  EXPECT_TRUE(boost::math::isinf(stan::math::squared_distance(v2, v1)));
}
