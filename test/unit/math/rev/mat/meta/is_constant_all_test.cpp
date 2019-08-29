#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename... Ts>
void expect_not_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_FALSE(temp);
}

using stan::length;

typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> var_t1;
typedef std::vector<var_t1> var_t2;
typedef std::vector<var_t2> var_t3;

typedef Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> var_u1;
typedef std::vector<var_u1> var_u2;
typedef std::vector<var_u2> var_u3;

typedef Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> var_v1;
typedef std::vector<var_v1> var_v2;
typedef std::vector<var_v2> var_v3;

TEST(MetaTraits, isConstantStruct) {
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
