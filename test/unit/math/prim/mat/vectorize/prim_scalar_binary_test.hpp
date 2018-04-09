#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_SCALAR_BINARY_TEST_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_SCALAR_BINARY_TEST_HPP

#include <test/unit/math/prim/mat/vectorize/expect_binary_types.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_values.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class prim_scalar_binary_test : public ::testing::Test {};

TYPED_TEST_CASE_P(prim_scalar_binary_test);

TYPED_TEST_P(prim_scalar_binary_test, expect_int_types) {
  expect_int_binary_types<TypeParam>();
}

TYPED_TEST_P(prim_scalar_binary_test, expect_scalar_types) {
  expect_binary_types<TypeParam, int, double>();
  expect_binary_types<TypeParam, double, int>();
  expect_binary_types<TypeParam, double, double>();
}

TYPED_TEST_P(prim_scalar_binary_test, expect_values) {
  expect_prim_binary_values<TypeParam>();
}

TYPED_TEST_P(prim_scalar_binary_test, expect_errors) {
  expect_prim_binary_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(prim_scalar_binary_test, expect_int_types,
                           expect_scalar_types, expect_values, expect_errors);
#endif
