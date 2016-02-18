#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_scalar_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_std_vector_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_matrix_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_vector_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_row_vector_value.hpp>

template <typename F>
void expect_values() {
  using stan::math::fvar;
  using stan::math::var;
  using std::vector;

  expect_scalar_value<F, fvar<var> >();
  expect_scalar_value<F, fvar<fvar<var> > >();
  expect_std_vector_value<F, fvar<var> >();
  expect_std_vector_value<F, fvar<fvar<var> > >();
  expect_matrix_value<F, fvar<var> >();
  expect_matrix_value<F, fvar<fvar<var> > >();
  expect_vector_value<F, fvar<var> >();
  expect_vector_value<F, fvar<fvar<var> > >();
  expect_row_vector_value<F, fvar<var> >();
  expect_row_vector_value<F, fvar<fvar<var> > >();
}

#endif
