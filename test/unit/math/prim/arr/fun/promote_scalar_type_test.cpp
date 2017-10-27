#include <gtest/gtest.h>
#include <stan/math/prim/arr.hpp>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>
#include <vector>

TEST(MathFunctionsPromoteScalarType, StdVector) {
  using std::vector;
  expect_promote_type<vector<double>, double, vector<int> >();
  expect_promote_type<vector<vector<double> >, double, vector<vector<int> > >();
  expect_promote_type<vector<vector<double> >, double,
                      vector<vector<double> > >();
}
