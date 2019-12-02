#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MetaTraitsPrimScal, ScalarTypeScal) {
  test::expect_same_type<double, stan::scalar_type<double>::type>();
  test::expect_same_type<int, stan::scalar_type<int>::type>();
}
