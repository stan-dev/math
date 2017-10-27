#include <gtest/gtest.h>
#include <stan/math/prim/arr.hpp>
#include <vector>

TEST(MetaTraits, containsNonconstantStruct) {
  using stan::contains_nonconstant_struct;
  using std::vector;
  EXPECT_FALSE(contains_nonconstant_struct<vector<double> >::value);
  EXPECT_FALSE(contains_nonconstant_struct<vector<vector<double> > >::value);
  EXPECT_FALSE(
      contains_nonconstant_struct<vector<vector<vector<double> > > >::value);
}
