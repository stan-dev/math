#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(ErrorHandlingMatrix, checkSizeMatch) {
  int x = 3;
  size_t y = 4;

  EXPECT_FALSE(stan::math::is_size_match(x, y));
  EXPECT_FALSE(stan::math::is_size_match(y, x));
  EXPECT_FALSE(stan::math::is_size_match(x, y));
  EXPECT_FALSE(stan::math::is_size_match(y, x));

  x = 2;
  y = 2;
  EXPECT_TRUE(stan::math::is_size_match(x, y));
  EXPECT_TRUE(stan::math::is_size_match(y, x));
  EXPECT_TRUE(stan::math::is_size_match(x, y));
  EXPECT_TRUE(stan::math::is_size_match(y, x));
}
