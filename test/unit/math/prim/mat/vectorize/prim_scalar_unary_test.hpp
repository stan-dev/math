#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_SCALAR_UNARY_TEST_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_SCALAR_UNARY_TEST_HPP

#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_values.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class prim_scalar_unary_test : public ::testing::Test {};

TYPED_TEST_CASE_P(prim_scalar_unary_test);

TYPED_TEST_P(prim_scalar_unary_test, expect_int_types) {
  expect_int_types<TypeParam>();
}

TYPED_TEST_P(prim_scalar_unary_test, expect_scalar_types) {
  expect_types<TypeParam, double>();
}

TYPED_TEST_P(prim_scalar_unary_test, expect_values) {
  expect_prim_values<TypeParam>();
}

TYPED_TEST_P(prim_scalar_unary_test, expect_errors) {
  expect_prim_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(prim_scalar_unary_test, expect_int_types,
                           expect_scalar_types, expect_values, expect_errors);
#endif
