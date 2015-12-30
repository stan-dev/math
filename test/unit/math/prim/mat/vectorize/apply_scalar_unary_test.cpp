#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_values.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

struct foo_base_test {
  static double constrain(double x) {
    return x;  // exp(x) for (0,inf), inv_logit(x) for (0,1), etc.
  }

  static std::vector<double> illegal_inputs() {
    return std::vector<double>();
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


// this tests that the expect_values test works on a mock function
TEST(MathPrimMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
}

// this tests that the return types work
template <typename F>
void expect_scalar_unary_return_type() {
  using stan::test::expect_match_return_t;
  expect_match_return_t<double, int>();

  typedef typename 
    stan::math::apply_scalar_unary<stan::math::foo_fun,
                                   std::vector<int> >::return_t
    vec_double_return_t;
  vec_double_return_t f;
  f.push_back(3.7);
  EXPECT_FLOAT_EQ(3.7, f[0]);

  expect_match_return_t<std::vector<double>, std::vector<int> >();
}
