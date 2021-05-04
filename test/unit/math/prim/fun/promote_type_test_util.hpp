#ifndef TEST_UNIT_MATH_PRIM_FUN_PROMOTE_TYPE_TEST_UTIL_HPP
#define TEST_UNIT_MATH_PRIM_FUN_PROMOTE_TYPE_TEST_UTIL_HPP

#include <stan/math/prim.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <type_traits>

template <typename T, typename S>
void expect_type(S s) {
  typedef BOOST_TYPEOF_TPL(stan::math::promote_scalar<T>(s)) result_t;
  bool same = std::is_same<S, result_t>::value;
  EXPECT_TRUE(same);
}

// pass:  expect_same_type<double, double>()
// fail:  expect_same_type<int, double>()
template <typename T, typename S>
void expect_same_type() {
  EXPECT_TRUE((std::is_same<S, T>::value));
}

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
