#include <stan/math/rev/meta.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::partials_return_type;
using stan::math::var;

TEST_F(AgradRev, MetaTraitsRevScal_PartialsReturnTypeVar) {
  EXPECT_SAME_TYPE(double, partials_return_type<var>::type);
}

TEST_F(AgradRev, MetaTraitsRevScal_PartialsReturnTypeVarTenParams) {
  EXPECT_SAME_TYPE(double,
                   partials_return_type<double, var, double, int, double, float,
                                        float, float, var, int>::type);
}

TEST_F(AgradRev, MetaTraitsRevArr_partials_return_type) {
  using stan::partials_return_type;
  using stan::math::var;

  partials_return_type<double, stan::math::var,
                       std::vector<stan::math::var> >::type g(5.0);
  EXPECT_EQ(5.0, g);
}
