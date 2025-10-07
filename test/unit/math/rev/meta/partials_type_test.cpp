#include <stan/math/rev/meta.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, MetaTraitsRevScal_partials_type) {
  using stan::partials_type;
  using stan::math::var;

  stan::partials_type<var>::type f(2.0);
  EXPECT_EQ(2.0, f);
}
