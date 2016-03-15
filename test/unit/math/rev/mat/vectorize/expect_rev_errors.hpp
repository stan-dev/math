#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_ERRORS_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_ERRORS_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_scalar_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_std_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_matrix_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_row_vector_error.hpp>

template <typename F>
void expect_rev_errors() {
  using stan::math::var;
  expect_scalar_error<F, var>();
  expect_std_vector_error<F, var>();
  expect_matrix_error<F, var>();
  expect_vector_error<F, var>();
  expect_row_vector_error<F, var>();
}

#endif
