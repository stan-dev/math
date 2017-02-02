#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::check_positive_finite;

TEST(AgradRevErrorHandlingScalar,CheckPositiveFinite_Vector) {
  const char* function = "check_positive_finite";
  std::vector<var> x;
  
  x.clear();
  x.push_back (1.5);
  x.push_back (0.1);
  x.push_back (1);
  ASSERT_NO_THROW(check_positive_finite(function, "x", x)) 
    << "check_positive_finite should be true with finite x";

  x.clear();
  x.push_back(1);
  x.push_back(2);
  x.push_back(std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error) 
    << "check_positive_finite should throw exception on Inf";

  x.clear();
  x.push_back(-1);
  x.push_back(2);
  x.push_back(std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error) 
    << "check_positive_finite should throw exception on negative x";

  x.clear();
  x.push_back(0);
  x.push_back(2);
  x.push_back(std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error) 
    << "check_positive_finite should throw exception on x=0";

  x.clear();
  x.push_back(1);
  x.push_back(2);
  x.push_back(-std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on -Inf";
  
  x.clear();
  x.push_back(1);
  x.push_back(2);
  x.push_back(std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on NaN";
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckPositiveFiniteVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_positive_finite;

  int N = 5;
  const char* function = "check_positive_finite";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U,stack_size);
  EXPECT_THROW(check_positive_finite(function,"a",a),std::domain_error);
  EXPECT_NO_THROW(check_positive_finite(function,"a",a[2]));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U,stack_size_after_call);

  a[2] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function,"a",a),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(6U,stack_size_after_call);
  stan::math::recover_memory();
}
