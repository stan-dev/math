#include <stan/math/prim/arr.hpp>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, value_type) {
  using stan::math::value_type;
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
