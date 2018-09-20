#include <gtest/gtest.h>
#include <stan/math/rev/scal.hpp>

TEST(MetaTraits, containsNonconstantStruct) {
  using stan::contains_nonconstant_struct;

  EXPECT_TRUE(contains_nonconstant_struct<stan::math::var>::value);
}
