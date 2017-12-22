#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, rising_factorial_var_int) {
  using stan::math::digamma;
  int a(1);
  AVAR b(4.0);
  AVAR f = stan::math::rising_factorial(b, a);
  EXPECT_FLOAT_EQ(4, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ((digamma(5.0) - digamma(4.0)) * 4.0, g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::rising_factorial(4.0 + eps, 1)
                   - stan::math::rising_factorial(4.0 - eps, 1))
                      / (2 * eps),
                  g[1]);
}

TEST(AgradRev, rising_factorial_exceptions) {
  int a(1);
  AVAR b(-3.0);
  EXPECT_NO_THROW(stan::math::rising_factorial(b, a));
  EXPECT_THROW(stan::math::rising_factorial(b, -1), std::domain_error);
}

struct rising_factorial_fun {
  template <typename T>
  inline typename stan::return_type<T>::type operator()(const T& arg1,
                                                        int arg2) const {
    return rising_factorial(arg1, arg2);
  }
};

TEST(AgradRev, check_varis_on_stack) {
  int a(1);
  stan::math::var b(4.0);

  test::check_varis_on_stack(stan::math::rising_factorial(b.val(), a));
  test::check_varis_on_stack(stan::math::rising_factorial(b.val(), 4));
  test::check_varis_on_stack(stan::math::rising_factorial(4, a));
}
