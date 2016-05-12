#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits,containsFvar) {
  using stan::math::fvar;
  using stan::contains_fvar;
  EXPECT_TRUE((contains_fvar<fvar<double> >::value));
  EXPECT_TRUE((contains_fvar<double, fvar<double> >::value));
  EXPECT_TRUE((contains_fvar<fvar<fvar<double> > >::value));
}
