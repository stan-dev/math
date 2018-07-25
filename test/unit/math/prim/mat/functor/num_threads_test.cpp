#define STAN_THREADS

#include <stan/math/prim/mat/functor/map_rect_concurrent.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>


#include <stdlib.h>

void set_n_threads_var(const std::string& value) {
  std::string env_string = "STAN_NUM_THREADS=" + value;
  putenv(env_string.c_str());
}

TEST(num_threads, correct_values) {
  set_n_threads_var("10");
  EXPECT_EQ(stan::math::internal::get_num_threads(100), 10);
  EXPECT_EQ(stan::math::internal::get_num_threads(5), 5);

  set_n_threads_var("-1");
  EXPECT_TRUE(stan::math::internal::get_num_threads(5) >= 1);
}

TEST(num_threads, incorrect_values) {
  set_n_threads_var("abc");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5), 
      std::runtime_error, "is not numeric");

  set_n_threads_var("1c");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5), 
      std::runtime_error, "is not numeric");

  set_n_threads_var("-2");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5), 
      std::runtime_error, "must be positive or -1");

  set_n_threads_var("0");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(5), 
      std::runtime_error, "must be positive or -1");

}