#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_SCALAR_UNARY_TEST_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_MIX_SCALAR_UNARY_TEST_HPP

#include <stan/math/mix/mat.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_types.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_values.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_mix_errors.hpp>

template <typename T>
class mix_scalar_unary_test : public ::testing::Test {
};

TYPED_TEST_CASE_P(mix_scalar_unary_test);

TYPED_TEST_P(mix_scalar_unary_test, expect_scalar_types) {
  using stan::math::fvar;
  using stan::math::var;

  expect_types<TypeParam, fvar<var> >();
  expect_types<TypeParam, fvar<fvar<var> > >();
}

TYPED_TEST_P(mix_scalar_unary_test, expect_values) {
  expect_mix_values<TypeParam>();
}

TYPED_TEST_P(mix_scalar_unary_test, expect_errors) {
  expect_mix_errors<TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P(mix_scalar_unary_test,
                           expect_scalar_types,
                           expect_values,
                           expect_errors);
#endif
