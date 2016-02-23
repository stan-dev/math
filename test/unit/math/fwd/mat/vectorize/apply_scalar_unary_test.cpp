#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/scal/fun/exp.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_values.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_errors.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_base_test.hpp>
#include <gtest/gtest.h>

// this tests that the expect_values test works on a mock function
TEST(MathFwdMatVectorize, applyScalarUnaryMock) {
  using stan::math::fvar;
  expect_types<foo_base_test, fvar<double> >();
  expect_types<foo_base_test, fvar<fvar<double> > >();
  expect_values<foo_base_test>();
  expect_errors<foo_base_test>();
}
