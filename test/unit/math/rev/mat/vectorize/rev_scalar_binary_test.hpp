#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_SCALAR_BINARY_TEST_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_SCALAR_BINARY_TEST_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_binary.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_types.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_values.hpp>
//#include <test/unit/math/rev/mat/vectorize/expect_rev_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class rev_scalar_binary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(rev_scalar_binary_test);

TYPED_TEST_P(rev_scalar_binary_test, expect_scalar_types) {
  using stan::math::var;
  expect_binary_types<TypeParam, var, int>();
  expect_binary_types<TypeParam, int, var>();
  expect_binary_types<TypeParam, var, double>();
  expect_binary_types<TypeParam, double, var>();
  expect_binary_types<TypeParam, var, var>();
}

TYPED_TEST_P(rev_scalar_binary_test, expect_values) {
  expect_rev_binary_values<TypeParam>();
}

/*
TYPED_TEST_P(rev_scalar_binary_test, expect_errors) {
  expect_rev_errors<TypeParam>();
}
*/
REGISTER_TYPED_TEST_CASE_P(rev_scalar_binary_test,
                           expect_scalar_types,
                           expect_values);
/*
                           expect_errors);
*/
#endif
