#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsNonconstantStruct) {
  using stan::contains_nonconstant_struct;

  EXPECT_TRUE(contains_nonconstant_struct<stan::math::var>::value);
}
