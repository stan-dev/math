#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isVectorMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  x.resize(3, 3);
  EXPECT_FALSE(stan::math::is_mat_vector(x));

  x.resize(0, 0);
  EXPECT_FALSE(stan::math::is_mat_vector(x));

  x.resize(1, 5);
  EXPECT_TRUE(stan::math::is_mat_vector(x));

  x.resize(5, 1);
  EXPECT_TRUE(stan::math::is_mat_vector(x));
}

TEST(ErrorHandlingMatrix, isVectorMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  double nan = std::numeric_limits<double>::quiet_NaN();

  x.resize(3, 3);
  x << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_FALSE(stan::math::is_mat_vector(x));

  x.resize(0, 0);
  EXPECT_FALSE(stan::math::is_mat_vector(x));

  x.resize(1, 5);
  x << nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_mat_vector(x));

  x.resize(5, 1);
  x << nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_mat_vector(x));
}
