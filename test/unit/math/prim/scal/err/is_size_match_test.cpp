#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(ErrorHandlingMatrix, checkSizeMatch) {
  using stan::math::is_size_match;
  int x = 3;
  size_t y = 4;

  EXPECT_FALSE(is_size_match(x, y));
  EXPECT_FALSE(is_size_match(y, x));
  EXPECT_FALSE(is_size_match(x, y)) << "int, size_t";
  EXPECT_FALSE(is_size_match(y, x)) << "size_t, int";

  x = 2;
  y = 2;
  EXPECT_TRUE(is_size_match(x, y));
  EXPECT_TRUE(is_size_match(y, x));
  EXPECT_TRUE(is_size_match(x, y)) << "int, size_t";
  EXPECT_TRUE(is_size_match(y, x)) << "size_t, int";
}
