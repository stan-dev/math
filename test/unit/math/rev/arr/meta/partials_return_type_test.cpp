#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, partials_return_type) {
  using stan::math::var;
  using stan::partials_return_type;

  partials_return_type<double,stan::math::var,std::vector<stan::math::var> >::type g(5.0);
  EXPECT_EQ(5.0,g);
}
