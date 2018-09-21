#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>

TEST(MetaTraits, get) {
  using stan::get;

  EXPECT_FLOAT_EQ(2.0, get(2.0, 1));
}
