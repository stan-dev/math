
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <vector>










TEST(MathFunctionsPromoteScalarType, primitive) {
  using std::vector;
  expect_promote_type<double, double, double>();
  expect_promote_type<double, double, int>();
}





TEST(MathFunctionsPromoteScalarType_arr, StdVector) {
  using std::vector;
  expect_promote_type<vector<double>, double, vector<int> >();
  expect_promote_type<vector<vector<double> >, double, vector<vector<int> > >();
  expect_promote_type<vector<vector<double> >, double,
                      vector<vector<double> > >();
}





TEST(MathFunctionsPromoteScalar_mat, TypeMatrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  expect_promote_type<Matrix<double, Dynamic, Dynamic>, double,
                      Matrix<int, Dynamic, Dynamic> >();

  expect_promote_type<Matrix<double, Dynamic, Dynamic>, double,
                      Matrix<double, Dynamic, Dynamic> >();

  expect_promote_type<vector<Matrix<double, Dynamic, Dynamic> >, double,
                      vector<Matrix<int, Dynamic, Dynamic> > >();
}

TEST(MathFunctionsPromoteScalar_mat, TypeVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  expect_promote_type<Matrix<double, Dynamic, 1>, double,
                      Matrix<int, Dynamic, 1> >();

  expect_promote_type<Matrix<double, Dynamic, 1>, double,
                      Matrix<double, Dynamic, 1> >();

  expect_promote_type<vector<Matrix<double, Dynamic, 1> >, double,
                      vector<Matrix<int, Dynamic, 1> > >();
}

TEST(MathFunctionsPromoteScalar_mat, TypeRowVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  expect_promote_type<Matrix<double, 1, Dynamic>, double,
                      Matrix<int, 1, Dynamic> >();

  expect_promote_type<Matrix<double, 1, Dynamic>, double,
                      Matrix<double, 1, Dynamic> >();

  expect_promote_type<vector<Matrix<double, 1, Dynamic> >, double,
                      vector<Matrix<int, 1, Dynamic> > >();
}
