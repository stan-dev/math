
#include <stan/math/prim.hpp>
#include <stan/math/prim/meta/index_type.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>
#include <vector>










TEST(MathMeta, index_type) {
  using stan::math::index_type;
  using std::vector;

  expect_same_type<vector<double>::size_type,
                   index_type<vector<double> >::type>();

  expect_same_type<vector<double>::size_type,
                   index_type<const vector<double> >::type>();

  expect_same_type<vector<vector<int> >::size_type,
                   index_type<vector<vector<int> > >::type>();

  expect_same_type<vector<vector<int> >::size_type,
                   index_type<const vector<vector<int> > >::type>();
}




TEST(MathMeta_mat, index_type) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::index_type;

  expect_same_type<Matrix<double, Dynamic, Dynamic>::Index,
                   index_type<Matrix<double, Dynamic, Dynamic> >::type>();

  expect_same_type<Matrix<double, Dynamic, 1>::Index,
                   index_type<Matrix<double, Dynamic, 1> >::type>();

  expect_same_type<Matrix<double, 1, Dynamic>::Index,
                   index_type<Matrix<double, 1, Dynamic> >::type>();
}
