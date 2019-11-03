#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_string_convertible) {
  using stan::is_string_convertible;
  EXPECT_TRUE(is_string_convertible<const char*>::value);
  EXPECT_FALSE(is_string_convertible<unsigned int>::value);
  EXPECT_FALSE(is_string_convertible<double>::value);
}
