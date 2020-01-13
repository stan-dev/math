#include <stan/math/prim/meta.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, value_type_vector) {
  using stan::value_type;
  using std::vector;

  expect_same_type<vector<double>::value_type,
                   value_type<vector<double> >::type>();

  expect_same_type<vector<double>::value_type,
                   value_type<const vector<double> >::type>();

  expect_same_type<vector<vector<int> >::value_type,
                   value_type<vector<vector<int> > >::type>();

  expect_same_type<vector<vector<int> >::value_type,
                   value_type<const vector<vector<int> > >::type>();
}

TEST(MathMetaPrim, value_type_matrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::value_type;

  expect_same_type<Matrix<double, Dynamic, Dynamic>::Scalar,
                   value_type<Matrix<double, Dynamic, Dynamic> >::type>();

  expect_same_type<Matrix<double, Dynamic, 1>::Scalar,
                   value_type<Matrix<double, Dynamic, 1> >::type>();

  expect_same_type<Matrix<double, 1, Dynamic>::Scalar,
                   value_type<Matrix<double, 1, Dynamic> >::type>();

  expect_same_type<
      Matrix<double, 1, Dynamic>,
      value_type<std::vector<Matrix<double, 1, Dynamic> > >::type>();
}
