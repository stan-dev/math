
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>




TEST(MetaTraits_scal, error_index) {
  EXPECT_EQ(1, static_cast<int>(stan::error_index::value));
}
