#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandling, validate_non_negative_index) {
  EXPECT_THROW(stan::math::validate_non_negative_index("unit_vector", "x", -5),
               std::invalid_argument);
  EXPECT_THROW(stan::math::validate_non_negative_index("unit_vector", "x", -1),
               std::invalid_argument);
  EXPECT_NO_THROW(
      stan::math::validate_non_negative_index("unit_vector", "x", 0));
}
