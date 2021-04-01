#include <stan/math/prim/core.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#ifdef TBB_INTERFACE_NEW

#include <tbb/global_control.h>
#include <tbb/task_arena.h>

TEST(intel_tbb_new_init, check_status) {
  set_n_threads(-1);
  tbb::task_arena& tbb_init = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(std::thread::hardware_concurrency(), tbb_init.max_concurrency());

  tbb::task_arena& tbb_reinit = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  tbb_init.terminate();
  EXPECT_FALSE(tbb_init.is_active());
}

#else

#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>

TEST(intel_tbb_init, check_status) {
  set_n_threads(-1);
  tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(std::thread::hardware_concurrency(),
            tbb::this_task_arena::max_concurrency());

  tbb::task_scheduler_init& tbb_reinit = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  tbb_init.terminate();
  EXPECT_FALSE(tbb_init.is_active());
}

#endif
