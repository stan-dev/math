#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_scalar_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_std_vector_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_matrix_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_vector_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_row_vector_value.hpp>

//Also tests derivatives
template <typename F>
void expect_values() {
  using stan::math::fvar;

  expect_scalar_value<F, fvar<double> >();
  expect_scalar_value<F, fvar<fvar<double> > >();
  expect_std_vector_value<F, fvar<double> >();
  expect_std_vector_value<F, fvar<fvar<double> > >();
  expect_matrix_value<F, fvar<double> >();
  expect_matrix_value<F, fvar<fvar<double> > >();
  expect_vector_value<F, fvar<double> >();
  expect_vector_value<F, fvar<fvar<double> > >();
  expect_row_vector_value<F, fvar<double> >();
  expect_row_vector_value<F, fvar<fvar<double> > >();
}

#endif
