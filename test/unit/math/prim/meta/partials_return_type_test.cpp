#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, PartialsReturnTypeDouble) {
  EXPECT_SAME_TYPE(double, stan::partials_return_type<double>::type);
}

TEST(MathMetaPrim, PartialsReturnTypeFloat) {
  EXPECT_SAME_TYPE(double, stan::partials_return_type<float>::type);
}

TEST(MathMetaPrim, PartialsReturnTypeInt) {
  EXPECT_SAME_TYPE(double, stan::partials_return_type<int>::type);
}

TEST(MathMetaPrim, PartialsReturnTypeScalarTenParams) {
  EXPECT_SAME_TYPE(
      double,
      stan::partials_return_type<double, int, double, float, float, double,
                                 float, int, double, double>::type);
}
