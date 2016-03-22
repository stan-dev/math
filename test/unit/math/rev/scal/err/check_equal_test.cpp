#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_equal;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar,CheckEqual) {
  const char* function = "check_equal";
  var x = 0.0;
  var eq = 0.0;
 
  EXPECT_TRUE(check_equal(function, "x", x, eq))
    << "check_equal should be true with x = eq";
  
  x = -1.0;
  EXPECT_THROW(check_equal(function, "x", x, eq),
               std::domain_error)
    << "check_equal should throw an exception with x < eq";

  x = eq;
  EXPECT_NO_THROW(check_equal(function, "x", x, eq))
    << "check_equal should not throw an exception with x == eq";

  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x, eq), 
               std::domain_error)
    << "check_equal should be false with x == Inf and eq = 0.0";

  x = 10.0;
  eq = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x, eq),
               std::domain_error)
    << "check_equal should throw an exception with x == 10.0 and eq == Inf";

  x = std::numeric_limits<double>::infinity();
  eq = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_equal(function, "x", x, eq))
    << "check_equal should not throw an exception with x == Inf and eq == Inf";
  stan::math::recover_memory();
}


TEST(AgradRevErrorHandlingScalar, CheckEqualVarCheckUnivariate) {

  const char* function = "check_equal";
  var a(5.0);
  var b(4.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(2U,stack_size);
  EXPECT_THROW(check_equal(function,"a",a,b),std::domain_error);

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(2U,stack_size_after_call);

  b = 5.0;
  EXPECT_TRUE(check_equal(function,"a",a,b));
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(3U,stack_size_after_call);

  b = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function,"a",a,b),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(4U,stack_size_after_call);
  stan::math::recover_memory();
}

