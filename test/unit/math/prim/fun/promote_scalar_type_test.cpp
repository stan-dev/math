#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctionsPromoteScalarType, primitive) {
  using std::vector;
  expect_promote_type<double, double, double>();
  expect_promote_type<double, double, int>();
}

TEST(MathFunctionsPromoteScalarType, StdVector) {
  using std::vector;
  expect_promote_type<vector<double>, double, vector<int> >();
  expect_promote_type<vector<vector<double> >, double, vector<vector<int> > >();
  expect_promote_type<vector<vector<double> >, double,
                      vector<vector<double> > >();
}

TEST(MathFunctionsPromoteScalar, TypeArray) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using std::vector;
  expect_promote_type<Array<double, Dynamic, Dynamic>, double,
                      Array<int, Dynamic, Dynamic> >();

  expect_promote_type<Array<double, Dynamic, Dynamic>, double,
                      Array<double, Dynamic, Dynamic> >();

  expect_promote_type<vector<Array<double, Dynamic, Dynamic> >, double,
                      vector<Array<int, Dynamic, Dynamic> > >();
}

TEST(MathFunctionsPromoteScalar, TypeMatrix) {
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

TEST(MathFunctionsPromoteScalar, TypeVector) {
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

TEST(MathFunctionsPromoteScalar, TypeRowVector) {
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

TEST(MathFunctionsPromoteScalar, TypeArrayExpression) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using std::vector;
  Array<int, Dynamic, Dynamic> a, b;
  using T_expr = decltype(a + b);
  expect_promote_type<Array<double, Dynamic, Dynamic>, double, T_expr>();
  expect_promote_type<Array<double, Dynamic, Dynamic>, double, T_expr>();
}

TEST(MathFunctionsPromoteScalar, TypeMatrixExpression) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  Matrix<int, Dynamic, Dynamic> a, b;
  using T_expr = decltype(a + b);
  expect_promote_type<Matrix<double, Dynamic, Dynamic>, double, T_expr>();
  expect_promote_type<Matrix<double, Dynamic, Dynamic>, double, T_expr>();
}
