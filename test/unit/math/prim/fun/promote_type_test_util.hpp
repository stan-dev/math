#ifndef TEST_UNIT_MATH_PRIM_FUN_PROMOTE_TYPE_TEST_UTIL_HPP
#define TEST_UNIT_MATH_PRIM_FUN_PROMOTE_TYPE_TEST_UTIL_HPP

#include <stan/math/prim.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <type_traits>

// E is expected value of promote_scalar_type<T, S>::type
template <typename E, typename T, typename S>
void expect_promote_type() {
  using stan::math::promote_scalar_type;
  using stan::math::test::type_name;
  using promoted_S = typename promote_scalar_type<T, S>::type;
  EXPECT_TRUE((std::is_same<E, promoted_S>::value))
      << "Expected Result: " << type_name<E>() << " from promoting "
      << type_name<S>() << " to " << type_name<T>() << " and instead got "
      << type_name<promoted_S>();
}

#endif
