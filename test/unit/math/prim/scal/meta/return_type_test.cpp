#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::return_type;

TEST(MetaTraitsPrimScal, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeScalarTenParams) {
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}
