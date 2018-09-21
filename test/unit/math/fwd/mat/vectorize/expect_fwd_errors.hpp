#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_ERRORS_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_ERRORS_HPP

#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_scalar_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_std_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_matrix_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_row_vector_error.hpp>

template <typename F>
void expect_fwd_errors() {
  using stan::math::fvar;
  expect_scalar_error<F, fvar<double> >();
  expect_scalar_error<F, fvar<fvar<double> > >();
  expect_std_vector_error<F, fvar<double> >();
  expect_std_vector_error<F, fvar<fvar<double> > >();
  expect_matrix_error<F, fvar<double> >();
  expect_matrix_error<F, fvar<fvar<double> > >();
  expect_vector_error<F, fvar<double> >();
  expect_vector_error<F, fvar<fvar<double> > >();
  expect_row_vector_error<F, fvar<double> >();
  expect_row_vector_error<F, fvar<fvar<double> > >();
}

#endif
