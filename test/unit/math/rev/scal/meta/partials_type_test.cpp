#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_type) {
  using stan::partials_type;
  using stan::math::var;

  stan::partials_type<var>::type f(2.0);
  EXPECT_EQ(2.0, f);
}
