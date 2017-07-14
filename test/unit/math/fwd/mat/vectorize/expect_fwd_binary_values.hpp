#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_VALUES_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_value.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_value.hpp>
#include <Eigen/Dense>

//Also will test derivatives
template <typename F>
void expect_fwd_binary_values() {
  using stan::math::fvar;
  expect_fwd_binary_scalar_value<F, fvar<double> >();
  expect_fwd_binary_std_vector_value<F, fvar<double> >();
  Eigen::MatrixXd model_matrix;
  Eigen::VectorXd model_vector;
  Eigen::RowVectorXd model_row_vector;
  expect_fwd_binary_matrix_value<F, fvar<double> >(model_matrix);
  expect_fwd_binary_matrix_value<F, fvar<double> >(model_vector);
  expect_fwd_binary_matrix_value<F, fvar<double> >(model_row_vector);
  expect_fwd_binary_scalar_value<F, fvar<fvar<double> > >();
  expect_fwd_binary_std_vector_value<F, fvar<fvar<double> > >();
  expect_fwd_binary_matrix_value<F, fvar<fvar<double> > >(model_matrix);
  expect_fwd_binary_matrix_value<F, fvar<fvar<double> > >(model_vector);
  expect_fwd_binary_matrix_value<F, fvar<fvar<double> > >(model_row_vector);
}

#endif
