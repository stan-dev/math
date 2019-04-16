#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingMatrix, isStdVectorIndexMatrix) {
  std::vector<double> y;
  y.push_back(5);
  y.push_back(5);
  y.push_back(5);
  y.push_back(5);
  size_t i;

  i = 2;
  y.resize(3);
  EXPECT_TRUE(stan::math::is_std_vector_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_std_vector_index(y, i));

  y.resize(2);
  EXPECT_FALSE(stan::math::is_std_vector_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_std_vector_index(y, i));
}

TEST(ErrorHandlingMatrix, isStdVectorIndexMatrix_nan) {
  std::vector<double> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  y.push_back(nan);
  y.push_back(nan);
  y.push_back(nan);
  y.push_back(nan);
  size_t i;

  i = 2;
  y.resize(3);
  EXPECT_TRUE(stan::math::is_std_vector_index(y, i));

  i = 3;
  EXPECT_TRUE(stan::math::is_std_vector_index(y, i));

  y.resize(2);
  EXPECT_FALSE(stan::math::is_std_vector_index(y, i));

  i = 0;
  EXPECT_FALSE(stan::math::is_std_vector_index(y, i));
}
