#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/arr/meta/index_type.hpp>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>
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
