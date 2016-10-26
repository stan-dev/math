#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,log_falling_factorial_var_double) {
  double a(2);
  AVAR b(4.0);
  AVAR f = stan::math::log_falling_factorial(b,a);
  EXPECT_FLOAT_EQ(std::log(12.0), f.val());

  AVEC x = createAVEC(a,b);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ((boost::math::digamma(5) - boost::math::digamma(3)),
                  g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(4.0 + eps, 2.0)
                  - stan::math::log_falling_factorial(4.0 - eps, 2.0))
                  / (2 * eps), g[1]);
}

TEST(AgradRev, log_falling_factorial_exceptions) {
  double a(1);
  AVAR b(-3.0);
  EXPECT_THROW(stan::math::log_falling_factorial(b,a), std::domain_error);
  EXPECT_THROW(stan::math::log_falling_factorial(b,b), std::domain_error);
}

TEST(AgradRev, log_falling_factorial_double_var) {
  double a(5);
  AVAR b(4.0);
  AVAR f = stan::math::log_falling_factorial(a,b);
  EXPECT_FLOAT_EQ(std::log(120.0), f.val());
  AVEC x = createAVEC(a,b);
  VEC g;  
  f.grad(x,g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ(boost::math::digamma(2), g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(5.0, 4.0 + eps)
                  - stan::math::log_falling_factorial(5.0, 4.0 - eps))
                  / (2 * eps), g[1]);
}

TEST(AgradRev, log_falling_factorial_var_var) {
  AVAR a(6.0);
  AVAR b(4.0);
  AVAR f = stan::math::log_falling_factorial(a,b);
  EXPECT_FLOAT_EQ(std::log(360.0), f.val());
  AVEC x = createAVEC(a,b);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ((boost::math::digamma(7) - boost::math::digamma(3)),
                  g[0]);
  EXPECT_FLOAT_EQ(boost::math::digamma(3), g[1]);
  
  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(6.0 + eps, 4.0)
                  - stan::math::log_falling_factorial(6.0 - eps, 4.0))
                  / (2 * eps), g[0]);
  EXPECT_FLOAT_EQ((stan::math::log_falling_factorial(6.0, 4.0 + eps)
                  - stan::math::log_falling_factorial(6.0, 4.0 - eps))
                  / (2 * eps), g[1]);
}

struct log_falling_factorial_fun {
  template <typename T0, typename T1>
  inline 
  typename stan::return_type<T0,T1>::type
  operator()(const T0& arg1,
             const T1& arg2) const {
    return log_falling_factorial(arg1,arg2);
  }
};

TEST(AgradRev, log_falling_factorial_nan) {
  log_falling_factorial_fun log_falling_factorial_;
  test_nan(log_falling_factorial_,4.0,4.0,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a(2.0);
  AVAR b(4.0);
  test::check_varis_on_stack(stan::math::log_falling_factorial(b, a));
  test::check_varis_on_stack(stan::math::log_falling_factorial(b, 2.0));
  test::check_varis_on_stack(stan::math::log_falling_factorial(4.0, a));
}
