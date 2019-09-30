#include <stan/math/prim/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/mat/functor/utils_threads.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>

TEST(intel_tbb_init, check_status) {
  const bool tbb_init = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init);

  EXPECT_EQ(tbb::task_scheduler_init::default_num_threads(),
            tbb::this_task_arena::max_concurrency());
}

TEST(intel_tbb_init, incorrect_env_values) {
  set_n_threads("abc");
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(true), std::invalid_argument,
                   "positive number or -1");

  set_n_threads("1c");
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(true), std::invalid_argument,
                   "positive number or -1");

  set_n_threads("-2");
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(true), std::invalid_argument,
                   "must be positive or -1");

  set_n_threads("0");
  EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(true), std::invalid_argument,
                   "must be positive or -1");
}
