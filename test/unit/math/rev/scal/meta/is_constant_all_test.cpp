#include <stan/math/rev/meta.hpp>
#include <gtest/gtest.h>

template <typename... Ts>
void expect_not_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_FALSE(temp);
}

TEST(MetaTraitsRevScal, isConstantStruct) {
  expect_not_const<stan::math::var>();
  expect_not_const<stan::math::var, double>();
  expect_not_const<stan::math::var, double, stan::math::var>();
  expect_not_const<stan::math::var, double, double>();
  expect_not_const<stan::math::var, stan::math::var, stan::math::var,
                   stan::math::var, stan::math::var, double, double>();
}

TEST(MetaTraitsRevArr, isConstantStruct) {
  using std::vector;

  expect_not_const<vector<stan::math::var> >();
  expect_not_const<vector<vector<stan::math::var> > >();
  expect_not_const<vector<vector<vector<stan::math::var> > > >();
  expect_not_const<vector<stan::math::var>, vector<double> >();
  expect_not_const<vector<stan::math::var>, vector<double>,
                   vector<stan::math::var> >();
  expect_not_const<vector<stan::math::var>,
                   vector<vector<vector<stan::math::var> > >,
                   vector<stan::math::var> >();
  expect_not_const<vector<stan::math::var>, vector<vector<vector<double> > >,
                   vector<stan::math::var> >();
}
