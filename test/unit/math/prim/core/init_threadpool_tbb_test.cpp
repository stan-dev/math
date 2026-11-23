#include <stan/math/prim/core.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

TEST(intel_tbb_new_init, check_status) {
  scoped_stan_num_threads restore_environment;
  set_n_threads(-1);
  auto& arena = stan::math::init_threadpool_tbb();
  ASSERT_TRUE(arena.is_active());
  const int original_threads = arena.max_concurrency();
#ifndef STAN_THREADS
  EXPECT_EQ(1, original_threads);
#endif

  // Other tests may have initialized the singleton before this test.
  set_n_threads(original_threads + 1);
  auto& reinit = stan::math::init_threadpool_tbb();
  EXPECT_EQ(&arena, &reinit);
  EXPECT_TRUE(reinit.is_active());
  EXPECT_EQ(original_threads, reinit.max_concurrency());
  auto& explicit_reinit = stan::math::init_threadpool_tbb(original_threads + 2);
  EXPECT_EQ(&arena, &explicit_reinit);
  EXPECT_EQ(original_threads, explicit_reinit.max_concurrency());
}

#ifdef STAN_THREADS
TEST(intel_tbb_new_init, validates_after_initialization) {
  scoped_stan_num_threads restore_environment;
  stan::math::init_threadpool_tbb(2);
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(-2), std::invalid_argument,
                   "positive or -1");
  set_n_threads("invalid");
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(), std::invalid_argument,
                   "positive number or -1");
}
#endif
