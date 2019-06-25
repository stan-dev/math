#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  using std::vector;
  EXPECT_TRUE(is_constant_all<vector<double>>::value);
  EXPECT_TRUE(is_constant_all<vector<vector<double>>>::value);
  EXPECT_TRUE(is_constant_all<vector<vector<vector<double>>>>::value);
  bool temp = is_constant_all<vector<double>, vector<double>,
                                 vector<double>>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_all<vector<double>, vector<vector<double>>,
                            vector<double>, vector<vector<double>>,
                            vector<vector<vector<double>>>>::value;
  EXPECT_TRUE(temp);
}
