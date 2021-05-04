#include <stan/math/rev/meta.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, partials_type) {
  using stan::partials_type;
  using stan::math::var;

  stan::partials_type<var>::type f(2.0);
  EXPECT_EQ(2.0, f);
}
