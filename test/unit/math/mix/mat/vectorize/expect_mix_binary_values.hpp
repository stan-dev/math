#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_VALUES_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_scalar_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_std_vector_value.hpp>
/*
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_matrix_value.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_vector_value.hpp>
*/
#include <Eigen/Dense>

//Also will test derivatives
template <typename F>
void expect_mix_binary_values() {
  using stan::math::var;
  using stan::math::fvar;
  expect_mix_binary_scalar_value<F, fvar<var> >();
  expect_mix_binary_std_vector_value<F, fvar<var> >();
/*
  expect_mix_binary_matrix_value<F, fvar<double> >();
  expect_mix_binary_vector_value<F, fvar<double>, Eigen::VectorXd>();
  expect_mix_binary_vector_value<F, fvar<double>, Eigen::RowVectorXd>();
  expect_mix_binary_scalar_value<F, fvar<fvar<double> > >();
  expect_mix_binary_std_vector_value<F, fvar<fvar<double> > >();
  expect_mix_binary_matrix_value<F, fvar<fvar<double> > >();
  expect_mix_binary_vector_value<F, fvar<fvar<double> >, Eigen::VectorXd>();
  expect_mix_binary_vector_value<F, fvar<fvar<double> >, 
  Eigen::RowVectorXd>();
*/
}

#endif
