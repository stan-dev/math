#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,rising_factorial_var_double) {
  using boost::math::digamma;
  double a(1);
  AVAR b(4.0);
  AVAR f = stan::math::rising_factorial(b,a);
  EXPECT_FLOAT_EQ(4,f.val());

  AVEC x = createAVEC(a,b);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ((digamma(5.0) - digamma(4.0)) * 4.0,g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(4.0 + eps, 1.0)
                  - stan::math::rising_factorial(4.0 - eps, 1.0))
                  / (2 * eps), g[1]);
} 

TEST(AgradRev, rising_factorial_exceptions) {
  double a(1);
  AVAR b(-3.0);
  EXPECT_THROW(stan::math::rising_factorial(b,a), std::domain_error);
  EXPECT_THROW(stan::math::rising_factorial(b,b), std::domain_error);
}

TEST(AgradRev, rising_factorial_double_var) {
  using boost::math::digamma;
  double a(5);
  AVAR b(4.0);
  AVAR f = stan::math::rising_factorial(a,b);
  EXPECT_FLOAT_EQ(5*6*7*8, f.val());
  AVEC x = createAVEC(a,b);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ(digamma(9) * 5*6*7*8, g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(5.0, 4.0 + eps)
                  - stan::math::rising_factorial(5.0, 4.0 - eps))
                  / (2 * eps), g[1]);
}

TEST(AgradRev, rising_factorial_var_var) {
  using boost::math::digamma;
  AVAR c(4.0);
  AVAR b(4.0);
  AVAR f = stan::math::rising_factorial(b,c);
  EXPECT_FLOAT_EQ(4*5*6*7, f.val());
  AVEC x = createAVEC(b,c);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(4.0*5.0*6.0*7.0 * (digamma(8.0) - digamma(4.0)), g[0]);
  EXPECT_FLOAT_EQ(4.0*5.0*6.0*7.0 * digamma(8), g[1]);
  
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(4.0 + eps, 4.0)
                  - stan::math::rising_factorial(4.0 - eps, 4.0))
                  / (2 * eps), g[0]);
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(4.0, 4.0 + eps)
                  - stan::math::rising_factorial(4.0, 4.0 - eps))
                  / (2 * eps), g[1]);
}

struct rising_factorial_fun {
  template <typename T0, typename T1>
  inline 
  typename stan::return_type<T0,T1>::type
  operator()(const T0& arg1,
             const T1& arg2) const {
    return rising_factorial(arg1,arg2);
  }
};

TEST(AgradRev, rising_factorial_nan) {
  rising_factorial_fun rising_factorial_;
  test_nan(rising_factorial_,3.0,5.0,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::var a(1);
  stan::math::var b(4.0);
  
  test::check_varis_on_stack(stan::math::rising_factorial(b, a));
  test::check_varis_on_stack(stan::math::rising_factorial(b, value_of(a)));
  test::check_varis_on_stack(stan::math::rising_factorial(value_of(b), a));
}
