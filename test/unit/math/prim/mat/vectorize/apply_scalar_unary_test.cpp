#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_values.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <cmath>
#include <iostream>

struct foo_base_test {
  typedef stan::math::foo_fun fun_t;

  static double constrain(double x) {
    return x;  // exp(x) for (0,inf), inv_logit(x) for (0,1), etc.
  }
  
  static bool has_illegal() {
    return false;
  }

  static double transform_illegal(double x) {
    return 0;
  }

  template <typename R, typename T>
  static R apply(const T& x) {
    using stan::math::foo;
    return foo(x);
  }

  static double apply_base(double x) {
    return apply<double>(x);
  }
};


TEST(MathPrimMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
}


