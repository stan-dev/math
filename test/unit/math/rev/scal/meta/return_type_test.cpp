#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::var;
using stan::return_type;

TEST(MetaTraits, ReturnTypeVar) {
  test::expect_same_type<var, return_type<var>::type>();
}

TEST(MetaTraits, ReturnTypeVarTenParams) {
  test::expect_same_type<var,
                         return_type<double, var, double, int, double, float,
                                     float, float, var, int>::type>();
}
