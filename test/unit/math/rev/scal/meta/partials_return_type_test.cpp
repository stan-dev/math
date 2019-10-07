#include <stan/math/rev/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::var;
using stan::partials_return_type;

TEST(MetaTraits, PartialsReturnTypeVar) {
  test::expect_same_type<double, partials_return_type<var>::type>();
}

TEST(MetaTraits, PartialsReturnTypeVarTenParams) {
  test::expect_same_type<
      double, partials_return_type<double, var, double, int, double, float,
                                   float, float, var, int>::type>();
}
