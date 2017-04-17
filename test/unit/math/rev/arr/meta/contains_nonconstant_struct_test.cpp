#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsNonconstantStruct) {
  using stan::contains_nonconstant_struct;
  using std::vector;

  EXPECT_TRUE(contains_nonconstant_struct<stan::math::var>::value);
  EXPECT_TRUE(contains_nonconstant_struct<vector<stan::math::var> >::value);
  EXPECT_TRUE(contains_nonconstant_struct<vector<vector<stan::math::var> > >::value);
  EXPECT_TRUE(contains_nonconstant_struct<vector<vector<vector<stan::math::var> > > >::value);
}
