#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_SCALAR_UNARY_TEST_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_FWD_SCALAR_UNARY_TEST_HPP

#include <stan/math/fwd/mat.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_values.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class fwd_scalar_unary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(fwd_scalar_unary_test);

TYPED_TEST_P(fwd_scalar_unary_test, expect_scalar_types) {
  using stan::math::fvar;

  expect_types<TypeParam, fvar<double> >();
  expect_types<TypeParam, fvar<fvar<double> > >();
}

TYPED_TEST_P(fwd_scalar_unary_test, expect_values) {
  expect_fwd_values<TypeParam>();
}

TYPED_TEST_P(fwd_scalar_unary_test, expect_errors) {
  expect_fwd_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(fwd_scalar_unary_test,
                           expect_scalar_types,
                           expect_values,
                           expect_errors);
#endif
