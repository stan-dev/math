#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::return_type;

TEST(MetaTraits, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraits, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraits, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraits, ReturnTypeScalarTenParams) {
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}
