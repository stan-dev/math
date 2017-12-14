#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_VALUES_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_VALUES_HPP

#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_scalar_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_std_vector_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_matrix_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_vector_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_row_vector_value.hpp>

// Also tests derivatives
template <typename F>
void expect_fwd_values() {
  using stan::math::fvar;

  expect_fwd_scalar_value<F, fvar<double> >();
  expect_fwd_scalar_value<F, fvar<fvar<double> > >();
  expect_fwd_std_vector_value<F, fvar<double> >();
  expect_fwd_std_vector_value<F, fvar<fvar<double> > >();
  expect_fwd_matrix_value<F, fvar<double> >();
  expect_fwd_matrix_value<F, fvar<fvar<double> > >();
  expect_fwd_vector_value<F, fvar<double> >();
  expect_fwd_vector_value<F, fvar<fvar<double> > >();
  expect_fwd_row_vector_value<F, fvar<double> >();
  expect_fwd_row_vector_value<F, fvar<fvar<double> > >();
}

#endif
