#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(ErrorHandlingMatrix, checkSizeMatch) {
  using stan::math::is_size_match;
  int x;
  size_t y;

  x = 3;
  y = 4;
  EXPECT_FALSE(is_size_match(x, y));
  EXPECT_FALSE(is_size_match(x, y));


  EXPECT_FALSE(is_size_match(x, y));
  EXPECT_FALSE(is_size_match(x, y));

  x = 2;
  y = 2;
  EXPECT_TRUE(is_size_match(x, y));
  EXPECT_TRUE(is_size_match(x, y));
  EXPECT_TRUE(is_size_match(x, y));
  EXPECT_TRUE(is_size_match(x, y));
}
