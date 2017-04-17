#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::check_positive_finite;

TEST(AgradRevErrorHandlingScalar,CheckPositiveFinite) {
  const char* function = "check_positive_finite";
  var x = 1;
 
  EXPECT_NO_THROW(check_positive_finite(function, "x", x))
    << "check_positive_finite should be true with finite x: " << x;
  x = -1;
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on x= " << x;
  x = 0;
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on x= " << x;
  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on Inf: " << x;
  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error) 
    << "check_positive_finite should throw exception on -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on NaN: " << x;
  stan::math::recover_memory();
}


TEST(AgradRevErrorHandlingScalar, CheckPositiveFiniteVarCheckUnivariate) {
  using stan::math::var;
  using stan::math::check_positive_finite;

  const char* function = "check_positive_finite";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(1U,stack_size);
  EXPECT_NO_THROW(check_positive_finite(function,"a",a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(1U,stack_size_after_call);

  stan::math::recover_memory();
}
