#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMetaFwd, value_type) {
  using stan::math::child_type;
  using stan::math::fvar;

  EXPECT_SAME_TYPE(double, child_type<fvar<double> >::type);
  EXPECT_SAME_TYPE(fvar<double>, child_type<fvar<fvar<double> > >::type);
}
