


#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
TEST(MetaTraits, error_index) {
  EXPECT_EQ(1, static_cast<int>(stan::error_index::value));
}
