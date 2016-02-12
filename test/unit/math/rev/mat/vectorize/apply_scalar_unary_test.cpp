#include <test/unit/math/prim/mat/vectorize/foo_base_test.hpp>
#include <test/unit/math/rev/mat/vectorize/rev_expect_values.hpp>
#include <test/unit/math/rev/mat/vectorize/rev_expect_errors.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <stan/math/rev/core/var.hpp>

// this tests that the expect_values test works on a mock function
TEST(MathRevMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
  if (foo_base_test::illegal_inputs().size() > 0)
    expect_errors<foo_base_test>();
}

// this tests that the return types work
template <typename F>
void expect_scalar_unary_return_type() {
   using stan::math::var;
   using stan::math::foo;
   using stan::test::expect_match_return_t;
   var three_var = 3;
   var exp_3_v = foo(three_var);
   EXPECT_FLOAT_EQ(std::exp(3.0), exp_3_v.val());
 
   expect_match_return_t<var, var>();
   expect_match_return_t<std::vector<var>, std::vector<var> >();
}
