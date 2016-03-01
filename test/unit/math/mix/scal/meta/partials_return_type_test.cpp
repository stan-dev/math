#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_return_type) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::partials_return_type;

  partials_return_type<double,fvar<fvar<var> > >::type c(3.0,2.0);
  EXPECT_EQ(3.0,c.val_.val());
  EXPECT_EQ(2.0,c.d_.val());

  partials_return_type<double,double,var,fvar<fvar<var> > >::type d(3.0,2.0);
  EXPECT_EQ(3.0,d.val_.val());
  EXPECT_EQ(2.0,d.d_.val());
}
