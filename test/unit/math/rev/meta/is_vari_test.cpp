#include <stan/math/rev/meta/is_vari.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_vari) {
  using stan::is_vari;
  EXPECT_TRUE(is_vari<stan::math::vari_value<double>>::value);
}
