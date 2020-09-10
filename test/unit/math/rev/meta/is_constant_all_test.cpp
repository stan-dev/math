#include <stan/math/rev/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

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

using var_t1 = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
using var_t2 = std::vector<var_t1>;
using var_t3 = std::vector<var_t2>;

using var_u1 = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
using var_u2 = std::vector<var_u1>;
using var_u3 = std::vector<var_u2>;

using var_v1 = Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic>;
using var_v2 = std::vector<var_v1>;
using var_v3 = std::vector<var_v2>;

TEST(MetaTraitsRevMat, isConstantStruct) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;

  expect_not_const<var_t1>();
  expect_not_const<var_t2>();
  expect_not_const<var_t3>();
  expect_not_const<var_u1>();
  expect_not_const<var_u2>();
  expect_not_const<var_u3>();
  expect_not_const<var_v1>();
  expect_not_const<var_v2>();
  expect_not_const<var_v3>();
  expect_not_const<var_t1, var_t2, var_t3, double>();
  expect_not_const<var_t1, var_u2, var_v3, double>();
  expect_not_const<var_t1, var_t2, var_t3, var_u1, var_u2, var_u3, var_v3,
                   double>();
}
