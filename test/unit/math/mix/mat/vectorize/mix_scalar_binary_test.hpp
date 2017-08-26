#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_SCALAR_BINARY_TEST_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_SCALAR_BINARY_TEST_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_types.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_values.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_binary_errors.hpp>
#include <gtest/gtest.h>

template <typename T>
class mix_scalar_binary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(mix_scalar_binary_test);

TYPED_TEST_P(mix_scalar_binary_test, expect_scalar_types) {
  using stan::math::var;
  using stan::math::fvar;
  expect_binary_types<TypeParam, fvar<var>, int>();
  expect_binary_types<TypeParam, int, fvar<var> >();
  expect_binary_types<TypeParam, fvar<var>, double>();
  expect_binary_types<TypeParam, double, fvar<var> >();
  expect_binary_types<TypeParam, fvar<var>, fvar<var> >();
  expect_binary_types<TypeParam, fvar<fvar<var> >, int>();
  expect_binary_types<TypeParam, int, fvar<fvar<var> > >();
  expect_binary_types<TypeParam, fvar<fvar<var> >, double>();
  expect_binary_types<TypeParam, double, fvar<fvar<var> > >();
  expect_binary_types<TypeParam, fvar<fvar<var> >, 
  fvar<fvar<var> > >();
}

TYPED_TEST_P(mix_scalar_binary_test, expect_values) {
  expect_mix_binary_values<TypeParam>();
}

TYPED_TEST_P(mix_scalar_binary_test, expect_errors) {
  expect_mix_binary_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(mix_scalar_binary_test,
                           expect_scalar_types,
                           expect_values,
                           expect_errors);
#endif
