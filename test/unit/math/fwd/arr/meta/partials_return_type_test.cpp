#include <stan/math/fwd/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, partials_return_type) {
  using stan::math::fvar;
  using stan::partials_return_type;

  partials_return_type<double, fvar<double>, std::vector<fvar<double> > >::type
      a(5.0);
  EXPECT_EQ(5.0, a);
}
