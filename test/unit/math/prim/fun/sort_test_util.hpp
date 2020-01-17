#ifndef TEST_UNIT_MATH_PRIM_FUN_SORT_TEST_UTIL_HPP
#define TEST_UNIT_MATH_PRIM_FUN_SORT_TEST_UTIL_HPP

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

template <typename C>
void test_sort_asc_throws() {
  using stan::math::sort_asc;

  C xs0;
  EXPECT_NO_THROW(sort_asc(xs0));

  C xs1(1);
  xs1[0] = 1;
  EXPECT_NO_THROW(sort_asc(xs1));
  xs1[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(sort_asc(xs1), std::domain_error);

  C xs2(2);
  xs2[0] = 1;
  xs2[1] = 2;
  EXPECT_NO_THROW(sort_asc(xs2));
  xs2[0] = 1;
  xs2[1] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(sort_asc(xs2), std::domain_error);
  xs2[0] = std::numeric_limits<double>::quiet_NaN();
  xs2[1] = 1;
  EXPECT_THROW(sort_asc(xs2), std::domain_error);
}

template <typename C>
void test_sort_desc_throws() {
  using stan::math::sort_desc;

  C xs0;
  EXPECT_NO_THROW(sort_desc(xs0));

  C xs1(1);
  xs1[0] = 1;
  EXPECT_NO_THROW(sort_desc(xs1));
  xs1[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(sort_desc(xs1), std::domain_error);

  C xs2(2);
  xs2[0] = 1;
  xs2[1] = 2;
  EXPECT_NO_THROW(sort_desc(xs2));
  xs2[0] = 1;
  xs2[1] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(sort_desc(xs2), std::domain_error);
  xs2[0] = std::numeric_limits<double>::quiet_NaN();
  xs2[1] = 1;
  EXPECT_THROW(sort_desc(xs2), std::domain_error);
}

#endif
