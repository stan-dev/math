#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, unitspaced_array) {
  using stan::math::unitspaced_array;
  EXPECT_STD_VECTOR_FLOAT_EQ(unitspaced_array(1, 5),
                             std::vector<int>({1, 2, 3, 4, 5}));
  EXPECT_STD_VECTOR_FLOAT_EQ(unitspaced_array(-2, 0),
                             std::vector<int>({-2, -1, 0}));
  EXPECT_STD_VECTOR_FLOAT_EQ(unitspaced_array(-5, -5), std::vector<int>({-5}));
}

TEST(MathFunctions, unitspaced_array_throw) {
  using stan::math::unitspaced_array;

  EXPECT_THROW(unitspaced_array(0, -1), std::domain_error);
  EXPECT_THROW(unitspaced_array(-5, -20), std::domain_error);
  EXPECT_THROW(unitspaced_array(50, 49), std::domain_error);
}
