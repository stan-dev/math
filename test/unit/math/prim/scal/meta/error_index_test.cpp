#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, error_index) {
  EXPECT_EQ(1, int(stan::error_index::value));
}
