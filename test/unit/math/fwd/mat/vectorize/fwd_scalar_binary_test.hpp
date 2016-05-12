#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_FWD_SCALAR_BINARY_TEST_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_FWD_SCALAR_BINARY_TEST_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_binary.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_types.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_values.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class fwd_scalar_binary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(fwd_scalar_binary_test);

TYPED_TEST_P(fwd_scalar_binary_test, expect_scalar_types) {
  using stan::math::fvar;
  expect_binary_types<TypeParam, fvar<double>, int>();
  expect_binary_types<TypeParam, int, fvar<double> >();
  expect_binary_types<TypeParam, fvar<double>, double>();
  expect_binary_types<TypeParam, double, fvar<double> >();
  expect_binary_types<TypeParam, fvar<double>, fvar<double> >();
  expect_binary_types<TypeParam, fvar<fvar<double> >, int>();
  expect_binary_types<TypeParam, int, fvar<fvar<double> > >();
  expect_binary_types<TypeParam, fvar<fvar<double> >, double>();
  expect_binary_types<TypeParam, double, fvar<fvar<double> > >();
  expect_binary_types<TypeParam, fvar<fvar<double> >, 
  fvar<fvar<double> > >();
}

TYPED_TEST_P(fwd_scalar_binary_test, expect_values) {
  expect_fwd_binary_values<TypeParam>();
}

TYPED_TEST_P(fwd_scalar_binary_test, expect_errors) {
  expect_fwd_binary_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(fwd_scalar_binary_test,
                           expect_scalar_types,
                           expect_values,
                           expect_errors);
#endif
