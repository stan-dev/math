#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, falling_factorial_var_int) {
  using stan::math::digamma;
  int a(2);
  AVAR b(4.0);
  AVAR f = stan::math::falling_factorial(b, a);
  EXPECT_FLOAT_EQ(12, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ((digamma(5) - digamma(3)) * 12.0, g[1]);

  double eps = 1e-6;
  EXPECT_FLOAT_EQ((stan::math::falling_factorial(4.0 + eps, 2)
                   - stan::math::falling_factorial(4.0 - eps, 2))
                      / (2 * eps),
                  g[1]);
}

TEST(AgradRev, falling_factorial_exceptions) {
  AVAR a = -1;
  int b(-3);
  EXPECT_THROW(stan::math::falling_factorial(a, b), std::domain_error);
  EXPECT_NO_THROW(stan::math::falling_factorial(b, 1));
}

struct falling_factorial_fun {
  template <typename T>
  inline typename stan::return_type<T>::type operator()(const T& arg1,
                                                        int arg2) const {
    return falling_factorial(arg1, arg2);
  }
};

TEST(AgradRev, check_varis_on_stack) {
  int a(2);
  AVAR b(4.0);
  test::check_varis_on_stack(stan::math::falling_factorial(b.val(), a));
  test::check_varis_on_stack(stan::math::falling_factorial(b.val(), 2));
  test::check_varis_on_stack(stan::math::falling_factorial(4, a));
}
