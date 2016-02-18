#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/scal/fun/exp.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_values.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_errors.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_base_test.hpp>

TEST(MathFwdMatVectorize, applyScalarUnaryMock) {
  expect_values<foo_base_test>();
  expect_errors<foo_base_test>();
}

