#include <stan/math/rev/meta.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_vari) {
  using stan::is_vari;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE(is_vari<stan::math::vari>::value);
  EXPECT_TRUE((is_vari<stan::math::vari_value<float>>::value));
  EXPECT_TRUE((is_vari<stan::math::vari_value<long double>>::value));
  EXPECT_FALSE(is_vari<stan::math::var>::value);
  EXPECT_FALSE((is_vari<double>::value));
  EXPECT_FALSE((is_vari<stan::math::var_value<float>>::value));
}
