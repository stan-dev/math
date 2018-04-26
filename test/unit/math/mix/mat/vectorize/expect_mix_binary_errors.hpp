#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_ERRORS_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_ERRORS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_std_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_matrix_error.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>

/**
 * Tests that the function defined statically in the template class
 * throws errors when supplied with illegal values.  The illegal
 * values are themselves defined statically in the template class.
 *
 * Tests integer, scalar, standard vector, Eigen matrix, Eigen vector,
 * and Eigen row vecor cases, and nested standard vectors of Eigen
 * types.
 *
 * @tparam F Test class used to define test case.
 */
template <typename F>
void expect_mix_binary_errors() {
  using stan::math::fvar;
  using stan::math::var;
  expect_binary_scalar_error<F, fvar<var> >();
  expect_binary_scalar_error<F, fvar<fvar<var> > >();
  expect_binary_std_vector_error<F, fvar<var> >();
  expect_binary_std_vector_error<F, fvar<fvar<var> > >();
  Eigen::MatrixXd model_matrix;
  Eigen::VectorXd model_vector;
  Eigen::RowVectorXd model_row_vector;
  expect_binary_matrix_error<F, fvar<var> >(model_matrix);
  expect_binary_matrix_error<F, fvar<fvar<var> > >(model_matrix);
  expect_binary_matrix_error<F, fvar<var> >(model_vector);
  expect_binary_matrix_error<F, fvar<fvar<var> > >(model_vector);
  expect_binary_matrix_error<F, fvar<var> >(model_row_vector);
  expect_binary_matrix_error<F, fvar<fvar<var> > >(model_row_vector);
}

#endif
