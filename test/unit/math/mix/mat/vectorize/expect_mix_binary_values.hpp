#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_VALUES_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_scalar_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_std_vector_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_matrix_value.hpp>
#include <Eigen/Dense>

// Also will test derivatives
template <typename F>
void expect_mix_binary_values() {
  using stan::math::fvar;
  using stan::math::var;
  expect_mix_binary_scalar_value<F, fvar<var> >();
  expect_mix_binary_scalar_value<F, fvar<fvar<var> > >();
  expect_mix_binary_std_vector_value<F, fvar<var> >();
  expect_mix_binary_std_vector_value<F, fvar<fvar<var> > >();

  Eigen::MatrixXd model_matrix;
  Eigen::VectorXd model_vector;
  Eigen::RowVectorXd model_row_vector;
  expect_mix_binary_matrix_value<F, fvar<var> >(model_matrix);
  expect_mix_binary_matrix_value<F, fvar<fvar<var> > >(model_matrix);
  expect_mix_binary_matrix_value<F, fvar<var> >(model_vector);
  expect_mix_binary_matrix_value<F, fvar<fvar<var> > >(model_vector);
  expect_mix_binary_matrix_value<F, fvar<var> >(model_row_vector);
  expect_mix_binary_matrix_value<F, fvar<fvar<var> > >(model_row_vector);
}

#endif
