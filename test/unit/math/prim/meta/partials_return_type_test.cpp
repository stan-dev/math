
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>






using stan::partials_return_type;

TEST(MetaTraits_scal, PartialsReturnTypeDouble) {
  test::expect_same_type<double, partials_return_type<double>::type>();
}

TEST(MetaTraits_scal, PartialsReturnTypeFloat) {
  test::expect_same_type<double, partials_return_type<float>::type>();
}

TEST(MetaTraits_scal, PartialsReturnTypeInt) {
  test::expect_same_type<double, partials_return_type<int>::type>();
}

TEST(MetaTraits_scal, PartialsReturnTypeScalarTenParams) {
  test::expect_same_type<
      double, partials_return_type<double, int, double, float, float, double,
                                   float, int, double, double>::type>();
}
