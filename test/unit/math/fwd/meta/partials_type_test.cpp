#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(MathMetaFwd, partials_type) {
  using stan::partials_type;
  using stan::math::fvar;

  stan::partials_type<fvar<double> >::type a(2.0);
  EXPECT_EQ(2.0, a);
  stan::partials_type<fvar<double> >::type b(4.0);
  EXPECT_EQ(4.0, b);
  stan::partials_type<fvar<fvar<double> > >::type c(7.0, 1.0);
  EXPECT_EQ(7.0, c.val_);
  EXPECT_EQ(1.0, c.d_);
}
