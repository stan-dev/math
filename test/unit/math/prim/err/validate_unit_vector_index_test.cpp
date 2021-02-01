#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandling, validate_unit_vector_index) {
  EXPECT_THROW(stan::math::validate_unit_vector_index("unit_vector", "x", 0),
               std::invalid_argument);
  EXPECT_THROW(stan::math::validate_unit_vector_index("unit_vector", "x", 1),
               std::invalid_argument);
  EXPECT_NO_THROW(
      stan::math::validate_unit_vector_index("unit_vector", "x", 5));
}
