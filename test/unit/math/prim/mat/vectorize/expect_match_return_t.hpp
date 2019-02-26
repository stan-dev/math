#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_MATCH_RETURN_T_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_MATCH_RETURN_T_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <gtest/gtest.h>
#include <type_traits>

namespace stan {

namespace test {

/**
 * Tests that the return type of a vectorized function with a
 * specified argument type is a specified type.
 *
 * @tparam F Functor defining function to test.
 * @tparam T_result_expected  Expected result type.
 * @tparam T_arg Argument type.
 */
template <typename F, typename T_result_expected, typename T_arg>
void expect_match_return_t() {
  using stan::math::apply_scalar_unary;
  typedef typename apply_scalar_unary<F, T_arg>::return_t result_t;
  EXPECT_TRUE((std::is_same<T_result_expected, result_t>::value));
}

}  // namespace test
}  // namespace stan
#endif
