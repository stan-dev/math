#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_SCALAR_UNARY_TEST_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_SCALAR_UNARY_TEST_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_values.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class rev_scalar_unary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(rev_scalar_unary_test);

TYPED_TEST_P(rev_scalar_unary_test, expect_scalar_types) {
  expect_types<TypeParam, stan::math::var>();
}

TYPED_TEST_P(rev_scalar_unary_test, expect_values) {
  expect_rev_values<TypeParam>();
}

TYPED_TEST_P(rev_scalar_unary_test, expect_errors) {
  expect_rev_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(rev_scalar_unary_test,
                           expect_scalar_types,
                           expect_values,
                           expect_errors);
#endif
