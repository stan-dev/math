#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::partials_return_type;

TEST(MetaTraitsPrimScal, PartialsReturnTypeDouble) {
  test::expect_same_type<double, partials_return_type<double>::type>();
}

TEST(MetaTraitsPrimScal, PartialsReturnTypeFloat) {
  test::expect_same_type<double, partials_return_type<float>::type>();
}

TEST(MetaTraitsPrimScal, PartialsReturnTypeInt) {
  test::expect_same_type<double, partials_return_type<int>::type>();
}

TEST(MetaTraitsPrimScal, PartialsReturnTypeScalarTenParams) {
  test::expect_same_type<
      double, partials_return_type<double, int, double, float, float, double,
                                   float, int, double, double>::type>();
}
