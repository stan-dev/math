#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename... Ts>
void expect_is_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_TRUE(temp);
}

TEST(MetaTraits, isConstantStruct) {
  using std::vector;
  expect_is_const<vector<double>>();
  expect_is_const<vector<vector<double>>>();
  expect_is_const<vector<vector<vector<double>>>>();
  expect_is_const<vector<double>, vector<double>, vector<double>>();
  expect_is_const<vector<double>, vector<vector<double>>, vector<double>,
                  vector<vector<double>>, vector<vector<vector<double>>>>();
}
