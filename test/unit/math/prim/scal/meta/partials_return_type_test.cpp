#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/scal.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::var;
using stan::partials_return_type;

TEST(MetaTraits, PartialsReturnTypeDouble) {
  test::expect_same_type<double, partials_return_type<double>::type>();
}

TEST(MetaTraits, PartialsReturnTypeFloat) {
  test::expect_same_type<double, partials_return_type<float>::type>();
}

TEST(MetaTraits, PartialsReturnTypeInt) {
  test::expect_same_type<double, partials_return_type<int>::type>();
}

TEST(MetaTraits, PartialsReturnTypeScalarTenParams) {
  test::expect_same_type<
      double, partials_return_type<double, int, double, float, float, double,
                                   float, int, double, double>::type>();
}
