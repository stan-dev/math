#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, PartialsReturnTypeDouble) {
  test::expect_same_type<double, stan::partials_return_type<double>::type>();
}

TEST(MathMetaPrim, PartialsReturnTypeFloat) {
  test::expect_same_type<double, stan::partials_return_type<float>::type>();
}

TEST(MathMetaPrim, PartialsReturnTypeInt) {
  test::expect_same_type<double, stan::partials_return_type<int>::type>();
}

TEST(MathMetaPrim, PartialsReturnTypeScalarTenParams) {
  test::expect_same_type<double, stan::partials_return_type<
                                     double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}
