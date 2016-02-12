#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_EXPECT_ERRORS_HPP

#include <test/unit/math/prim/mat/vectorize/expect_errors.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/rev/core.hpp>

template <typename F>
void expect_errors() {
  using stan::math::fvar;
  using stan::math::var;

  expect_scalar_error<F, fvar<var> >();
  expect_scalar_error<F, fvar<fvar<var> > >();
  expect_std_vectors_error<F, fvar<var> >();
  expect_std_vectors_error<F, fvar<fvar<var> > >();
  expect_matrix_error<F, fvar<var> >();
  expect_matrix_error<F, fvar<fvar<var> > >();
  expect_vector_error<F, fvar<var> >();
  expect_vector_error<F, fvar<fvar<var> > >();
  expect_row_vector_error<F, fvar<var> >();
  expect_row_vector_error<F, fvar<fvar<var> > >();
}

#endif
