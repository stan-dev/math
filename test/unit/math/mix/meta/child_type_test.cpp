#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMetaMix, value_type) {
  using stan::math::child_type;
  using stan::math::fvar;
  using stan::math::var;

  EXPECT_SAME_TYPE(fvar<var>, child_type<fvar<fvar<var> > >::type);
  EXPECT_SAME_TYPE(var, child_type<fvar<var> >::type);
}
