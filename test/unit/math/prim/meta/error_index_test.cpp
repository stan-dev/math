#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsPrimScal, error_index) {
  EXPECT_EQ(1, static_cast<int>(stan::error_index::value));
}
