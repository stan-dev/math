#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>
#include <gtest/gtest.h>

TEST(get_num_threads, correct_values) {
  set_n_threads("10");
  EXPECT_EQ(stan::math::internal::get_num_threads(), 10);

  set_n_threads("4");
  EXPECT_EQ(stan::math::internal::get_num_threads(), 4);

  set_n_threads("-1");
  EXPECT_EQ(stan::math::internal::get_num_threads(),
            std::thread::hardware_concurrency());
}

TEST(get_num_threads, incorrect_values) {
  set_n_threads("abc");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(),
                   std::invalid_argument, "positive number or -1");

  set_n_threads("1c");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(),
                   std::invalid_argument, "positive number or -1");

  set_n_threads("-2");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(),
                   std::invalid_argument, "must be positive or -1");

  set_n_threads("0");
  EXPECT_THROW_MSG(stan::math::internal::get_num_threads(),
                   std::invalid_argument, "must be positive or -1");
}
