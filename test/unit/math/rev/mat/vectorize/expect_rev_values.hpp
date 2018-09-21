#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_VALUES_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_VALUES_HPP

#include <test/unit/math/rev/mat/vectorize/expect_rev_scalar_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_std_vector_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_matrix_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_vector_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_row_vector_value.hpp>

// Also will test derivatives
template <typename F>
void expect_rev_values() {
  expect_rev_scalar_value<F>();
  expect_rev_std_vector_value<F>();
  expect_rev_matrix_value<F>();
  expect_rev_vector_value<F>();
  expect_rev_row_vector_value<F>();
}

#endif
