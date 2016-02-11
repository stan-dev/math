#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_values.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_errors.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <stan/math/fwd/core/fvar.hpp>

struct foo_base_test {

  template <typename R, typename T>
  static R apply(const T& x) {
    using stan::math::foo;
    return foo(x);
  }

  template <typename T>
  static stan::math::fvar<T> apply_base(stan::math::fvar<T> x) {
    return apply<stan::math::fvar<T> >(x);
  }

  static std::vector<double> valid_inputs() {
    using std::vector;

    vector<double> valid_inputs;
    valid_inputs.push_back(1.3);
    valid_inputs.push_back(-2.6);
    valid_inputs.push_back(0);
    valid_inputs.push_back(-0.2);

    return valid_inputs;
  }

  static std::vector<double> illegal_inputs() {
    using std::vector;

    vector<double> illegal_inputs(2, 10.6);
    illegal_inputs.push_back(25.7);
    illegal_inputs.push_back(100.25);

    return illegal_inputs;
  }
};


// this tests that the expect_values test works on a mock function
TEST(MathFwdMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
  if (foo_base_test::illegal_inputs().size() > 0)
    expect_errors<foo_base_test>();
}

// this tests that the return types work
template <typename F>
void expect_scalar_unary_return_type() {
   using stan::math::fvar;
   using stan::math::foo;
   using stan::test::expect_match_return_t;
   fvar<double> three_var = 3;
   fvar<double> exp_3_v = foo(three_var);
   EXPECT_FLOAT_EQ(std::exp(3.0), exp_3_v.val());
 
   expect_match_return_t<fvar<double>, fvar<double> >();
   expect_match_return_t<std::vector<fvar<double> >, 
                           std::vector<fvar<double> > >();
}
