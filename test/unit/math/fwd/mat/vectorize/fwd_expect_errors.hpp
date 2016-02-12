#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_EXPECT_ERRORS_HPP

#include <stdexcept>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stan/math/fwd/core.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_errors.hpp>

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_errors() {
  using stan::math::fvar;

  expect_scalar_error<F, fvar<double> >();
  expect_scalar_error<F, fvar<fvar<double> > >();
  expect_std_vectors_error<F, fvar<double> >();
  expect_std_vectors_error<F, fvar<fvar<double> > >();
  expect_matrix_error<F, fvar<double> >();
  expect_matrix_error<F, fvar<fvar<double> > >();
  expect_vector_error<F, fvar<double> >();
  expect_vector_error<F, fvar<fvar<double> > >();
  expect_row_vector_error<F, fvar<double> >();
  expect_row_vector_error<F, fvar<fvar<double> > >();
}

#endif
