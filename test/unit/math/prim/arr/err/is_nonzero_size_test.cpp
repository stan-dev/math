#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <vector>

TEST(ErrorHandlingMatrix, IsNonzeroSizeMatrix) {
  using stan::math::is_nonzero_size;

  std::vector<double> a;
  a.push_back(3.0);
  a.push_back(3.0);
  a.push_back(3.0);
  a.push_back(3.0);

  EXPECT_TRUE(stan::math::is_nonzero_size(a));

  a.resize(2);
  EXPECT_TRUE(stan::math::is_nonzero_size(a));

  a.resize(0);
  EXPECT_FALSE(stan::math::is_nonzero_size(a));
}

TEST(ErrorHandlingMatrix, IsNonzeroSizeMatrix_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> a;
  a.push_back(nan);
  a.push_back(nan);
  a.push_back(nan);
  a.push_back(nan);

  EXPECT_TRUE(stan::math::is_nonzero_size( a));

  a.resize(2);
  EXPECT_TRUE(stan::math::is_nonzero_size( a));

  a.resize(0);
  EXPECT_FALSE(stan::math::is_nonzero_size( a));
}
