#include <test/unit/math/prim/mat/vectorize/expect_values.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_errors.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_base_test.hpp>
#include <gtest/gtest.h>

// test that returned values are correct and that errors are thrown
// for illegal inputs for foo
TEST(MathPrimMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
  expect_errors<foo_base_test>();
}
