#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(MetaTraits, is_string_convertible) {
  using stan::is_string_convertible;
  EXPECT_TRUE(is_string_convertible<const char*>::value);
  EXPECT_TRUE(is_string_convertible<std::string>::value);
  EXPECT_FALSE(is_string_convertible<unsigned int>::value);
  EXPECT_FALSE(is_string_convertible<double>::value);
}
