#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  using std::vector;

  EXPECT_FALSE(is_constant_struct<vector<stan::math::var> >::value);
  EXPECT_FALSE(is_constant_struct<vector<vector<stan::math::var> > >::value);
  EXPECT_FALSE(
      is_constant_struct<vector<vector<vector<stan::math::var> > > >::value);
}
