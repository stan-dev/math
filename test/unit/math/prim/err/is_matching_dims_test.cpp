#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, isMatchingDimsMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  y.resize(3, 3);
  x.resize(3, 3);
  EXPECT_TRUE(stan::math::is_matching_dims(x, y));

  x.resize(0, 0);
  y.resize(0, 0);
  EXPECT_TRUE(stan::math::is_matching_dims(x, y));

  y.resize(1, 2);
  EXPECT_FALSE(stan::math::is_matching_dims(x, y));

  x.resize(2, 1);
  EXPECT_FALSE(stan::math::is_matching_dims(x, y));
}

TEST(ErrorHandlingMatrix, isMatchingDimsMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  x.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  x << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_TRUE(stan::math::is_matching_dims(x, y));

  x.resize(0, 0);
  y.resize(0, 0);
  EXPECT_TRUE(stan::math::is_matching_dims(x, y));

  y.resize(1, 2);
  y << nan, nan;
  EXPECT_FALSE(stan::math::is_matching_dims(x, y));

  x.resize(2, 1);
  x << nan, nan;
  EXPECT_FALSE(stan::math::is_matching_dims(x, y));
}

TEST(ErrorHandlingMatrix, isMatchingDims_compile_time_sizes) {
  using stan::math::is_matching_dims;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_dynamic;
  Eigen::Matrix<double, 2, 2> m_2x2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> vector(4);
  Eigen::Matrix<double, 1, Eigen::Dynamic> rowvector(4);

  m_dynamic.resize(2, 2);
  EXPECT_TRUE(is_matching_dims(m_dynamic, m_2x2));

  m_dynamic.resize(3, 3);
  EXPECT_FALSE(is_matching_dims(m_dynamic, m_2x2));

  m_dynamic.resize(4, 1);
  EXPECT_TRUE(is_matching_dims(m_dynamic, vector));

  m_dynamic.resize(3, 1);
  EXPECT_FALSE(is_matching_dims(m_dynamic, vector));

  m_dynamic.resize(1, 4);
  EXPECT_TRUE(is_matching_dims(m_dynamic, rowvector));

  m_dynamic.resize(1, 3);
  EXPECT_FALSE(is_matching_dims(m_dynamic, rowvector));
}
