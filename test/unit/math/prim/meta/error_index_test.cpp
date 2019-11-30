#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsPrimScal, error_index) {
  EXPECT_EQ(1, static_cast<int>(stan::error_index::value));
}
