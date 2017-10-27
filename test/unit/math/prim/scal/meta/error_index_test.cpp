#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>

TEST(MetaTraits, error_index) {
  EXPECT_EQ(1, static_cast<int>(stan::error_index::value));
}
