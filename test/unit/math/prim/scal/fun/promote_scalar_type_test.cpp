#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>

TEST(MathFunctionsPromoteScalarType, primitive) {
  using std::vector;
  expect_promote_type<double, double, double>();
  expect_promote_type<double, double, int>();
}
