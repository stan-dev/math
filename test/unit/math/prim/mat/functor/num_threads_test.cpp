#ifndef STAN_THREADS
#define STAN_THREADS
#endif

#include <stan/math/prim/mat/functor/map_rect_concurrent.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

#include <stdlib.h>

// Can't easily use std::string as putenv require non-const char*
void set_n_threads_var(const char* value) {
  static char env_string[256];
  snprintf(env_string, sizeof(env_string), "STAN_NUM_THREADS=%s", value);
  putenv(env_string);
}

TEST(num_threads, correct_values) {
  set_n_threads_var("10");
  EXPECT_EQ(stan::math::internal::get_num_threads(100), 10);

  set_n_threads_var("4");
  EXPECT_EQ(stan::math::internal::get_num_threads(3), 3);

  set_n_threads_var("-1");
  EXPECT_GE(stan::math::internal::get_num_threads(5), 1);
}

TEST(num_threads, incorrect_values) {
  set_n_threads_var("abc");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5),
                   std::invalid_argument, "positive number or -1");

  set_n_threads_var("1c");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5),
                   std::invalid_argument, "positive number or -1");

  set_n_threads_var("-2");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5),
                   std::invalid_argument, "must be positive or -1");

  set_n_threads_var("0");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5),
                   std::invalid_argument, "must be positive or -1");
}
