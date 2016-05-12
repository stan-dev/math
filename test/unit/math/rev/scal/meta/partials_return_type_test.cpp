#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_return_type) {
  using stan::math::var;
  using stan::partials_return_type;

  partials_return_type<double,stan::math::var>::type f(5.0);
  EXPECT_EQ(5.0,f);
}
