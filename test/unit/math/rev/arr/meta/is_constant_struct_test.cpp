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
  bool temp
      = is_constant_struct<vector<stan::math::var>, vector<double> >::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<vector<stan::math::var>, vector<double>,
                            vector<stan::math::var> >::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<vector<stan::math::var>,
                            vector<vector<vector<stan::math::var> > >,
                            vector<stan::math::var> >::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<vector<stan::math::var>,
                            vector<vector<vector<double> > >,
                            vector<stan::math::var> >::value;
  EXPECT_FALSE(temp);
}
