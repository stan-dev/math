#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  using stan::math::fvar;
  EXPECT_TRUE((contains_fvar<fvar<double> >::value));
  EXPECT_TRUE((contains_fvar<double, fvar<double> >::value));
  EXPECT_TRUE((contains_fvar<fvar<fvar<double> > >::value));
  bool temp = contains_fvar<double, double, double, double, double,
                            fvar<fvar<double> > >::value;
  EXPECT_TRUE(temp);
}
