#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  using std::vector;
  EXPECT_TRUE(is_constant_struct<vector<double> >::value);
  EXPECT_TRUE(is_constant_struct<vector<vector<double> > >::value);
  EXPECT_TRUE(is_constant_struct<vector<vector<vector<double> > > >::value);
}
