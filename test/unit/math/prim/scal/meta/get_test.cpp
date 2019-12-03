#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsPrimScal, get) {
  using stan::get;

  EXPECT_FLOAT_EQ(2.0, get(2.0, 1));
}
