#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VALUES_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VALUES_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_scalar_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_std_vector_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_matrix_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_vector_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_row_vector_value.hpp>

template <typename F>
void expect_values() {
  expect_scalar_value<F>();
  expect_std_vector_value<F>();
  expect_matrix_value<F>();
  expect_vector_value<F>();
  expect_row_vector_value<F>();
}

#endif
