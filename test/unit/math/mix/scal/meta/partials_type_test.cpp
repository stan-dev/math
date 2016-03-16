#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_type) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::partials_type;

  stan::partials_type<fvar<fvar<var> > >::type d(7.0,1.0);
  EXPECT_EQ(7.0,d.val_.val());
  EXPECT_EQ(1.0,d.d_.val());
  stan::partials_type<fvar<var> >::type e(2.0);
  EXPECT_EQ(2.0,e.val());
}
