#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsNonconstantStruct) {
  using stan::contains_nonconstant_struct;
  EXPECT_FALSE(contains_nonconstant_struct<int>::value);
  EXPECT_FALSE(contains_nonconstant_struct<double>::value);
  EXPECT_FALSE(contains_nonconstant_struct<float>::value);
  EXPECT_FALSE(contains_nonconstant_struct<int32_t>::value);
}
