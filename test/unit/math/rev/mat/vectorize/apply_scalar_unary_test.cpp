#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_base_test.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_values.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_errors.hpp>
#include <gtest/gtest.h>

/**
 * this tests that the expect_values test works on a mock function 
 */
TEST(MathRevMatVectorize, applyScalarUnaryMock) {
  expect_types<foo_base_test, stan::math::var>();
  expect_values<foo_base_test>();
  expect_errors<foo_base_test>();
}

