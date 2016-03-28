#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_return_type) {
  using stan::math::fvar;
  using stan::partials_return_type;

  partials_return_type<double,fvar<fvar<double> > >::type b(3.0,2.0);
  EXPECT_EQ(3.0,b.val_);
  EXPECT_EQ(2.0,b.d_);

  partials_return_type<double,double, fvar<fvar<double> >, fvar<fvar<double> > >::type e(3.0,2.0);
  EXPECT_EQ(3.0,e.val_);
  EXPECT_EQ(2.0,e.d_);

}
