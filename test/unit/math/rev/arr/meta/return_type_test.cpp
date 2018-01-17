#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

using stan::math::var;
using stan::return_type;
using std::vector;

TEST(MetaTraits, ReturnTypeVarArray) {
  test::expect_same_type<var, return_type<vector<var> >::type>();
  test::expect_same_type<var, return_type<vector<var>, double>::type>();
  test::expect_same_type<var, return_type<vector<var>, double>::type>();
}

TEST(MetaTraits, ReturnTypeDoubleArray) {
  test::expect_same_type<double, return_type<vector<double> >::type>();
  test::expect_same_type<double, return_type<vector<double>, double>::type>();
  test::expect_same_type<double, return_type<vector<double>, double>::type>();
}
