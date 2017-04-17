#include <stan/math/prim/arr.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, ScalarTypeArray) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<vector<double> >::type>();
  test::expect_same_type<int, scalar_type<vector<int> >::type>();
  test::expect_same_type<double, scalar_type<vector<vector<double> > >::type>();
}
