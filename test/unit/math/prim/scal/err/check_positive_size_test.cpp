#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>
#include <string>
#include <test/unit/util.hpp>

TEST(ErrorHandlingScalar, CheckPositiveSize) {
  using stan::math::check_positive_size;
  const std::string function = "function";
  const std::string name = "name";
  const std::string expr = "expr";
  std::string expected_msg;

  EXPECT_NO_THROW(check_positive_size(function, name, expr, 10));

  expected_msg =
      "name must have a positive size, but is 0; "
      "dimension size expression = expr";
  EXPECT_THROW_MSG(check_positive_size(function, name, expr, 0),
                   std::invalid_argument, expected_msg);

  expected_msg =
      "name must have a positive size, but is -1; "
      "dimension size expression = expr";
  EXPECT_THROW_MSG(check_positive_size(function, name, expr, -1),
                   std::invalid_argument, expected_msg);
}
