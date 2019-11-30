#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename... Ts>
void expect_is_const() {
  using stan::is_constant_all;
  bool temp = is_constant_all<Ts...>::value;
  EXPECT_TRUE(temp);
}
TEST(MetaTraitsPrimScal, isConstantStruct) {
  expect_is_const<>();
  expect_is_const<int>();
  expect_is_const<double>();
  expect_is_const<float>();
  expect_is_const<unsigned int>();
  expect_is_const<int32_t>();
  expect_is_const<int, int>();
  expect_is_const<double, double, double>();
  expect_is_const<float, float, float, float>();
  expect_is_const<int32_t, int32_t, int32_t, int32_t>();
}

TEST(MetaTraitsPrimArr, isConstantStruct) {
  using std::vector;
  expect_is_const<vector<double>>();
  expect_is_const<vector<vector<double>>>();
  expect_is_const<vector<vector<vector<double>>>>();
  expect_is_const<vector<double>, vector<double>, vector<double>>();
  expect_is_const<vector<double>, vector<vector<double>>, vector<double>,
                  vector<vector<double>>, vector<vector<vector<double>>>>();
}

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const_t1;
typedef std::vector<const_t1> const_t2;
typedef std::vector<const_t2> const_t3;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> const_u1;
typedef std::vector<const_u1> const_u2;
typedef std::vector<const_u2> const_u3;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> const_v1;
typedef std::vector<const_v1> const_v2;
typedef std::vector<const_v2> const_v3;

TEST(MetaTraitsPrimMat, isConstantStruct) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  expect_is_const<const_t1>();
  expect_is_const<const_t2>();
  expect_is_const<const_t3>();
  expect_is_const<const_u1>();
  expect_is_const<const_u2>();
  expect_is_const<const_u3>();
  expect_is_const<const_v1>();
  expect_is_const<const_v2>();
  expect_is_const<const_v3>();
  expect_is_const<const_t1, const_t2>();
  expect_is_const<const_t2, const_t3, const_u1>();
  expect_is_const<const_u1, const_v1, const_v2, const_t2>();
}
