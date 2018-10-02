#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_VALUES_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_VALUES_HPP

#include <stan/math/mix/mat.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_scalar_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_std_vector_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_matrix_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_vector_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_row_vector_value.hpp>

// Also will test derivatives
template <typename F>
void expect_mix_values() {
  using stan::math::fvar;
  using stan::math::var;
  using std::vector;

  expect_mix_scalar_value<F, fvar<var> >();
  expect_mix_scalar_value<F, fvar<fvar<var> > >();
  expect_mix_std_vector_value<F, fvar<var> >();
  expect_mix_std_vector_value<F, fvar<fvar<var> > >();
  expect_mix_matrix_value<F, fvar<var> >();
  expect_mix_matrix_value<F, fvar<fvar<var> > >();
  expect_mix_vector_value<F, fvar<var> >();
  expect_mix_vector_value<F, fvar<fvar<var> > >();
  expect_mix_row_vector_value<F, fvar<var> >();
  expect_mix_row_vector_value<F, fvar<fvar<var> > >();
}

#endif
